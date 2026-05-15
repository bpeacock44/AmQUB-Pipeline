#!/usr/bin/env python3
"""
Per-OTU Spearman correlation between
  - nematode-sample abundance (Set A)
  - soil fold-change after/before treatment (Set B) computed across SoilNumbers.

For each OTU we produce:
  - nem_s   = per-soil nematode abundance:
                within soil s, mean each --setA group separately,
                then mean those group-means
  - fc_s    = per-soil fold change:
                (mean(--setB_after samples in soil s) + c) /
                (mean(--setB_before samples in soil s) + c)

Then Spearman(nem_vec, fc_vec) across soils → one rho/p per OTU.
BH-FDR is applied across all OTUs.
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
    parser = argparse.ArgumentParser(description="Per-OTU TU-FC correlation across soils")
    parser.add_argument('--asv_file', required=True)
    parser.add_argument('--metadata_file', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--column', required=True, help='Group column in metadata')
    parser.add_argument('--treatment_col', required=True, help='Soil/treatment column in metadata')
    parser.add_argument('--setA', required=True, help='Comma-separated groups for Set A (nematodes)')
    parser.add_argument('--setB_before', required=True, help='Comma-separated groups for Set B (before)')
    parser.add_argument('--setB_after', required=True, help='Comma-separated groups for Set B (after)')
    parser.add_argument('--pseudocount', type=float, default=1.0,
                        help='Pseudocount added before computing fold change (default 1)')
    return parser.parse_args()


# --------------------------------------------------
# Logging
# --------------------------------------------------

def setup_logger(output_dir, column):
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f'corr_otu_v_foldincr_{column}_{current_time}.log')
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

def load_asv_table(asv_file):
    with open(asv_file, 'r') as f:
        first_line = f.readline().strip()
        skiprows = [0] if first_line == "# Constructed from biom file" else []
    return pd.read_csv(asv_file, sep='\t', index_col=0, skiprows=skiprows)

def load_metadata(metadata_file):
    return pd.read_csv(metadata_file, sep='\t', index_col=0)

# --------------------------------------------------
# Aggregation
# --------------------------------------------------

def aggregate(asv, meta, group_col, treatment_col):
    """Mean abundance per (Group, SoilNumber) cell. Returns DataFrame
    indexed by OTU with MultiIndex columns (group, soil)."""
    meta = meta.dropna(subset=[group_col])
    meta = meta.loc[meta.index.intersection(asv.columns)]
    out = {}
    for (group, soil), samples in meta.groupby([group_col, treatment_col]):
        sub = asv[samples.index]
        if sub.empty:
            continue
        out[(group, soil)] = sub.mean(axis=1)
    return pd.DataFrame(out)


def per_soil_average(agg_df, groups):
    """Select MultiIndex columns whose group is in `groups`, then mean
    across groups within each soil. Returns DataFrame indexed by OTU
    with one column per soil. Realises: mean per group, then mean of means."""
    sel = agg_df.loc[:, agg_df.columns.get_level_values(0).isin(groups)]
    if sel.empty:
        return pd.DataFrame()
    return sel.T.groupby(level=1, sort=False).mean().T


# --------------------------------------------------
# Per-OTU correlation across soils
# --------------------------------------------------

def per_otu_correlation(nem_by_soil, fc_by_soil, outdir, column_name):
    # Align on common soils and OTUs
    soils = nem_by_soil.columns.intersection(fc_by_soil.columns)
    if len(soils) == 0:
        raise ValueError("No SoilNumbers in common between Set A and Set B fold-change")
    nem_by_soil = nem_by_soil[soils]
    fc_by_soil = fc_by_soil[soils]

    otus = nem_by_soil.index.intersection(fc_by_soil.index)
    nem_by_soil = nem_by_soil.loc[otus]
    fc_by_soil = fc_by_soil.loc[otus]

    n_soils = len(soils)
    logging.info(f"Per-OTU Spearman across {n_soils} soils: {list(soils)}")
    if n_soils < 4:
        logging.warning(
            f"Only {n_soils} soils — Spearman with n<4 is unreliable "
            "(ties give NaN; otherwise rho is bounded to a few discrete values "
            "and p-values are uninformative)."
        )

    soil_str = ','.join(map(str, soils))
    rows = []
    for otu in otus:
        x = nem_by_soil.loc[otu].to_numpy(dtype=float)
        y = fc_by_soil.loc[otu].to_numpy(dtype=float)
        rho, p = spearmanr(x, y, nan_policy="omit")
        rows.append({
            'TU': otu,
            'Correlation Coefficient': rho,
            'P-Value': p,
            'Avg_Abun': float(np.mean(x)),
            'Fold_Increase': float(np.mean(y)),
            'TU_Counts_SetA': ','.join(f'{v:g}' for v in x),
            'TU_Counts_SetB_Ratio': ','.join(f'{v:g}' for v in y),
            'Sample_IDs': soil_str,
        })
    df = pd.DataFrame(rows)

    # BH-FDR across OTUs (skip NaN p-values)
    pvals = df['P-Value'].to_numpy()
    fdr = np.full_like(pvals, np.nan, dtype=float)
    mask = ~np.isnan(pvals)
    if mask.sum() > 0:
        _, fdr_corrected, _, _ = multipletests(pvals[mask], method='fdr_bh')
        fdr[mask] = fdr_corrected
    df.insert(3, 'FDR', fdr)

    # Write outputs
    os.makedirs(outdir, exist_ok=True)
    out_txt = os.path.join(outdir, f"otu_correlation_results_{column_name}.txt")
    out_xlsx = out_txt.replace('.txt', '.xlsx')
    df.to_csv(out_txt, sep='\t', index=False)
    df.to_excel(out_xlsx, index=False)
    return df


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    setup_logger(args.output_dir, args.column)

    logging.info("Loading data")
    asv = load_asv_table(args.asv_file)
    meta = load_metadata(args.metadata_file)

    agg = aggregate(asv, meta, args.column, args.treatment_col)
    if agg.empty:
        raise ValueError("No data after aggregation")

    setA_groups = [g.strip() for g in args.setA.split(',')]
    setB_before_groups = [g.strip() for g in args.setB_before.split(',')]
    setB_after_groups = [g.strip() for g in args.setB_after.split(',')]
    logging.info(f"Set A groups:        {setA_groups}")
    logging.info(f"Set B before groups: {setB_before_groups}")
    logging.info(f"Set B after groups:  {setB_after_groups}")

    nem_by_soil    = per_soil_average(agg, setA_groups)
    before_by_soil = per_soil_average(agg, setB_before_groups)
    after_by_soil  = per_soil_average(agg, setB_after_groups)

    if nem_by_soil.empty:
        raise ValueError(f"No samples for Set A groups: {setA_groups}")
    if before_by_soil.empty:
        raise ValueError(f"No samples for Set B (before) groups: {setB_before_groups}")
    if after_by_soil.empty:
        raise ValueError(f"No samples for Set B (after) groups: {setB_after_groups}")

    common_soils = before_by_soil.columns.intersection(after_by_soil.columns)
    if len(common_soils) == 0:
        raise ValueError("No SoilNumbers in common between Set B before and after")
    c = args.pseudocount
    fc_by_soil = (after_by_soil[common_soils] + c) / (before_by_soil[common_soils] + c)

    per_otu_correlation(nem_by_soil, fc_by_soil, args.output_dir, args.column)
    logging.info("Done")


if __name__ == "__main__":
    main()