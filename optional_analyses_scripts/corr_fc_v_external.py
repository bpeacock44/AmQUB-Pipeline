#!/usr/bin/env python3
"""
Per-OTU Spearman correlation between
  - per-soil values from a prior output of corr_otu_v_foldincr.py
    (by default, the per-soil fold change in `TU_Counts_SetB_Ratio`)
  - an unrelated per-soil value supplied as a 2-column TSV (soil, value).

For each OTU we:
  1. Read its per-soil value vector and matching soil IDs from the prior output.
  2. Look up the external value for each soil.
  3. Drop soils with NaN on either side.
  4. Spearman(otu_vec, external_vec) across the remaining soils.

BH-FDR is applied across all OTUs. Output mirrors corr_otu_v_foldincr.py.
"""

import argparse
import os
import logging
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


# --------------------------------------------------
# Arguments
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Per-OTU correlation between prior per-soil FC and an external per-soil value"
    )
    parser.add_argument('--prior_results', required=True,
                        help='TSV output from corr_otu_v_foldincr.py')
    parser.add_argument('--external_file', required=True,
                        help='Two-column TSV: first column = soil ID, second = value. '
                             'Header is auto-detected.')
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--value_column', default='TU_Counts_SetB_Ratio',
                        help='Which per-soil vector column to use from the prior results '
                             '(default: TU_Counts_SetB_Ratio = per-soil fold change). '
                             'Use TU_Counts_SetA to correlate nematode abundance instead.')
    parser.add_argument('--value_name', default='External',
                        help='Short label for the external value, used in output column names '
                             '(e.g. pH -> Avg_pH, pH_Values). Default: External.')
    return parser.parse_args()


# --------------------------------------------------
# Logging
# --------------------------------------------------

def setup_logger(output_dir, basename):
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f'{basename}_corr_to_external_{current_time}.log')
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger = logging.getLogger()
    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
               for h in logger.handlers):
        logger.addHandler(console)


# --------------------------------------------------
# IO
# --------------------------------------------------

def load_prior_results(path, value_column):
    df = pd.read_csv(path, sep='\t')
    for required in ('TU', value_column, 'Sample_IDs'):
        if required not in df.columns:
            raise ValueError(
                f"Required column '{required}' not in prior results. "
                f"Columns present: {list(df.columns)}"
            )
    return df


def load_external_values(path):
    """Two-column TSV: soil, value. Auto-detect header by sniffing row 1.
    Returns dict {soil_id_str: float_value}."""
    raw = pd.read_csv(path, sep='\t', header=None, dtype=str)
    if raw.shape[1] < 2:
        raise ValueError(
            f"External file must have at least 2 columns (soil, value); got {raw.shape[1]}"
        )
    # Sniff: try to coerce second cell of row 0 to float. If it fails, row 0 is a header.
    first_val = raw.iloc[0, 1]
    has_header = False
    try:
        float(first_val)
    except (ValueError, TypeError):
        has_header = True

    if has_header:
        soil_col, val_col = raw.iloc[0, 0], raw.iloc[0, 1]
        logging.info(f"External file: detected header (soil='{soil_col}', value='{val_col}')")
        data = raw.iloc[1:].copy()
    else:
        logging.info("External file: no header detected")
        data = raw

    soils = data.iloc[:, 0].astype(str).str.strip()
    values = pd.to_numeric(data.iloc[:, 1], errors='coerce')

    n_bad = int(values.isna().sum())
    if n_bad > 0:
        bad_rows = data.loc[values.isna()].iloc[:, :2].values.tolist()
        logging.warning(
            f"External file: {n_bad} row(s) had non-numeric values and were dropped: {bad_rows}"
        )

    mapping = dict(zip(soils, values))
    # Drop NaN entries from the map so per-OTU lookups can treat "missing soil" uniformly.
    mapping = {k: v for k, v in mapping.items() if pd.notna(v)}

    if len(mapping) == 0:
        raise ValueError("External file contained no usable (soil, value) pairs")
    logging.info(f"External file: {len(mapping)} soils with usable values")
    return mapping


# --------------------------------------------------
# Helpers
# --------------------------------------------------

def parse_csv_floats(s):
    """Parse a comma-separated string of numbers into a numpy float array.
    Non-numeric tokens (incl. 'nan', empty) become NaN."""
    if pd.isna(s) or str(s).strip() == '':
        return np.array([], dtype=float)
    out = []
    for tok in str(s).split(','):
        tok = tok.strip()
        if tok == '' or tok.lower() == 'nan':
            out.append(np.nan)
            continue
        try:
            out.append(float(tok))
        except ValueError:
            out.append(np.nan)
    return np.array(out, dtype=float)


def parse_csv_strs(s):
    if pd.isna(s):
        return []
    return [tok.strip() for tok in str(s).split(',')]


# --------------------------------------------------
# Per-OTU correlation
# --------------------------------------------------

def per_otu_correlation(prior_df, external_map, value_column, value_name,
                         outdir, basename):
    rows = []
    used_soils_global = set()
    n_skipped_no_overlap = 0
    n_low_n = 0

    for _, r in prior_df.iterrows():
        otu = r['TU']
        soils = parse_csv_strs(r['Sample_IDs'])
        vals = parse_csv_floats(r[value_column])

        if len(soils) != len(vals):
            logging.warning(
                f"OTU {otu}: Sample_IDs length ({len(soils)}) != "
                f"{value_column} length ({len(vals)}); skipping."
            )
            continue

        # Build paired vectors restricted to soils present in the external map
        x, y, used = [], [], []
        for soil, v in zip(soils, vals):
            ext = external_map.get(soil)
            if ext is None or pd.isna(v) or pd.isna(ext):
                continue
            x.append(v)
            y.append(ext)
            used.append(soil)

        if len(used) == 0:
            n_skipped_no_overlap += 1
            continue
        used_soils_global.update(used)

        if len(used) < 4:
            n_low_n += 1

        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        if len(used) >= 2:
            rho, p = spearmanr(x, y)
        else:
            rho, p = np.nan, np.nan

        rows.append({
            'TU': otu,
            'Correlation Coefficient': rho,
            'P-Value': p,
            f'Avg_{value_column}': float(np.mean(x)) if len(x) else np.nan,
            f'Avg_{value_name}': float(np.mean(y)) if len(y) else np.nan,
            'N_Soils_Used': len(used),
            f'{value_column}_Values': ','.join(f'{v:g}' for v in x),
            f'{value_name}_Values': ','.join(f'{v:g}' for v in y),
            'Sample_IDs': ','.join(used),
        })

    if len(rows) == 0:
        raise ValueError(
            "No OTUs produced a correlation. Check that soil IDs in the external file "
            "match those in the prior results' Sample_IDs column."
        )

    df = pd.DataFrame(rows)
    logging.info(
        f"Computed correlation for {len(df)} OTUs; "
        f"{n_skipped_no_overlap} skipped (no soil overlap with external file); "
        f"{n_low_n} had <4 paired soils."
    )
    logging.info(f"Soils used across at least one OTU: {sorted(used_soils_global)}")

    if n_low_n > 0:
        logging.warning(
            f"{n_low_n} OTU(s) had fewer than 4 paired soils. Spearman with n<4 is unreliable "
            "(ties give NaN; otherwise rho is bounded to a few discrete values "
            "and p-values are uninformative)."
        )

    # BH-FDR
    pvals = df['P-Value'].to_numpy()
    fdr = np.full_like(pvals, np.nan, dtype=float)
    mask = ~np.isnan(pvals)
    if mask.sum() > 0:
        _, fdr_corrected, _, _ = multipletests(pvals[mask], method='fdr_bh')
        fdr[mask] = fdr_corrected
    df.insert(3, 'FDR', fdr)

    os.makedirs(outdir, exist_ok=True)
    out_txt = os.path.join(outdir, f"{basename}_corr_to_external.txt")
    out_xlsx = os.path.join(outdir, f"{basename}_corr_to_external.xlsx")
    df.to_csv(out_txt, sep='\t', index=False)
    df.to_excel(out_xlsx, index=False)
    logging.info(f"Wrote {out_txt}")
    logging.info(f"Wrote {out_xlsx}")
    return df


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    basename = os.path.splitext(os.path.basename(args.prior_results))[0]
    setup_logger(args.output_dir, basename)

    logging.info(f"Prior results:  {args.prior_results}")
    logging.info(f"External file:  {args.external_file}")
    logging.info(f"Value column:   {args.value_column}")
    logging.info(f"Value name:     {args.value_name}")

    prior = load_prior_results(args.prior_results, args.value_column)
    external = load_external_values(args.external_file)

    per_otu_correlation(prior, external, args.value_column, args.value_name,
                         args.output_dir, basename)
    logging.info("Done")


if __name__ == "__main__":
    main()