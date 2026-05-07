#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
import datetime

"""
Unified ASV/OTU processing script
"""

# ---------------------------
# Logger
# ---------------------------
def setup_logger(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = os.path.join(output_dir, f'prism_maker_{current_time}.log')

    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


# ---------------------------
# Helpers
# ---------------------------
def safe_name(name):
    return str(name).replace(" ", "_").replace("/", "_")

def prepare_abundance(df, context=""):
    df_num = df.apply(pd.to_numeric, errors='coerce')
    num_na = df_num.isna().sum().sum()
    if num_na > 0:
        logging.warning(f"{num_na} values were coerced to NaN during numeric conversion ({context})")
    all_nan_cols = df_num.columns[df_num.isna().all()]
    if len(all_nan_cols) > 0:
        logging.warning(f"{len(all_nan_cols)} columns are entirely NaN after coercion ({context})")
    return df_num

def build_row(avg_abundance, taxa_list, treatment_name):
    row = {'treatment': treatment_name}
    values = avg_abundance['average_abundance']
    for t in taxa_list:
        row[t] = values.get(t, 0)
    row['other'] = values.drop(index=taxa_list, errors='ignore').sum()
    total = sum(v for k, v in row.items() if k != 'treatment')
    if total == 0:
        logging.warning(f"Total abundance is zero for treatment '{treatment_name}'")
    if total > 0:
        for k in row:
            if k != 'treatment':
                row[k] = (row[k] / total) * 100
    return row

def sort_columns_by_abundance(df):
    if df.empty:
        logging.warning("Attempted to sort columns on empty dataframe")
        return df
    value_cols = [c for c in df.columns if c not in ['treatment', 'other']]
    sorted_taxa = sorted(
        value_cols,
        key=lambda x: df.iloc[0][x],
        reverse=True
    )
    ordered_cols = ['treatment'] + sorted_taxa + ['other']
    return df[ordered_cols]


# ---------------------------
# Arguments
# ---------------------------
parser = argparse.ArgumentParser()

parser.add_argument('--file_path', required=True)
parser.add_argument('--map_file', required=True)
parser.add_argument('--treatment_column', required=True)
parser.add_argument('--filter_column', required=True)
parser.add_argument('--output_dir', required=True)
parser.add_argument('--num_taxa', type=int, default=3)

parser.add_argument(
    '--mode',
    choices=['average', 'per_treatment'],
    default='average'
)

args = parser.parse_args()

setup_logger(args.output_dir)

# ---------------------------
# Load ASV table
# ---------------------------
logging.info(f"Reading ASV table: {args.file_path}")

with open(args.file_path) as f:
    first_line = f.readline().strip()

skip_rows = 1 if first_line == "# Constructed from biom file" else 0

df = pd.read_csv(args.file_path, sep='\t', skiprows=skip_rows, comment='~')
df.rename(columns={df.columns[0]: "ASV_ID"}, inplace=True)

# Validate ASV structure
if 'taxonomy' not in df.columns and 'ASV_ID' not in df.columns:
    raise ValueError("ASV table must contain a 'taxonomy' column or a valid first column.")

if 'taxonomy' in df.columns:
    df['taxonomy'] = df['taxonomy'] + '_' + df['ASV_ID']
    df.drop(columns=['ASV_ID'], inplace=True)
else:
    df.rename(columns={'ASV_ID': 'taxonomy'}, inplace=True)

# ---------------------------
# Load mapping file
# ---------------------------
logging.info(f"Reading mapping file: {args.map_file}")

map_df = pd.read_csv(args.map_file, sep='\t', comment='~')
map_df.rename(columns={'#SampleID': 'sampleID'}, inplace=True)

# Validate mapping file columns
required_map_cols = {'sampleID', args.treatment_column, args.filter_column}
missing_cols = required_map_cols - set(map_df.columns)
if missing_cols:
    raise ValueError(f"Mapping file missing required columns: {missing_cols}")

# Keep only relevant columns
map_df = map_df[['sampleID', args.treatment_column, args.filter_column]]

# Filter samples
map_df = map_df[~pd.isna(map_df[args.filter_column])]
map_df = map_df[~pd.isna(map_df[args.treatment_column])]

treatments = map_df[args.treatment_column].unique()
logging.info(f"Number of treatments: {len(treatments)}")

# ---------------------------
# Prepare ASV table
# ---------------------------
df_t = df.set_index('taxonomy').T
df_t.reset_index(inplace=True)
df_t.rename(columns={'index': 'sampleID'}, inplace=True)

# Warn about mismatched samples
pre_merge_samples = set(map_df['sampleID'])
df_samples = set(df_t['sampleID'])

missing_in_asv = pre_merge_samples - df_samples
missing_in_map = df_samples - pre_merge_samples

if missing_in_asv:
    logging.warning(f"{len(missing_in_asv)} samples in mapping file not found in ASV table")

if missing_in_map:
    logging.warning(f"{len(missing_in_map)} samples in ASV table not found in mapping file")

# Merge
merged_df = pd.merge(map_df, df_t, on='sampleID')

if merged_df.empty:
    logging.warning("Merged dataframe is empty after join")

# ---------------------------
# Overall abundance average
# ---------------------------
overall_abundance = merged_df.drop(
    columns=['sampleID', args.treatment_column, args.filter_column]
)

overall_abundance = prepare_abundance(overall_abundance, context="overall")

overall_avg = (
    overall_abundance.mean(axis=0)
    .rename_axis('taxonomy')
    .reset_index(name='average_abundance')
    .set_index('taxonomy')
)

# ---------------------------
# Compute averages + top taxa
# ---------------------------
treatment_avg = {}
unique_top_taxa = set()

for treatment in treatments:
    logging.info(f"Processing treatment: {treatment}")
    subset = merged_df[merged_df[args.treatment_column] == treatment]
    if subset.empty:
        logging.warning(f"Treatment '{treatment}' has no samples after filtering")
        continue
    abundance = subset.drop(columns=['sampleID', args.treatment_column, args.filter_column])
    # unified numeric cleanup (same as overall)
    abundance = prepare_abundance(abundance, context=f"treatment={treatment}")
    avg = (
        abundance.mean(axis=0)
        .rename_axis('taxonomy')
        .reset_index(name='average_abundance')
        .sort_values(by='average_abundance', ascending=False)
        .set_index('taxonomy')
    )
    treatment_avg[treatment] = avg
    top_taxa = avg.head(args.num_taxa).index.tolist()
    unique_top_taxa.update(top_taxa)

final_top_taxa = sorted(unique_top_taxa)
logging.info(f"Total unique top taxa: {len(final_top_taxa)}")

# ---------------------------
# Average output
# ---------------------------
if args.mode == 'average':
    logging.info("Generating global average abundance PRISM file")
    logging.info(f"Output file will be called: average_top{args.num_taxa}_{args.treatment_column}_{args.filter_column}.tsv")
    row = build_row(overall_avg, final_top_taxa, "average")
    output_df = pd.DataFrame([row])
    output_df = sort_columns_by_abundance(output_df)
    filename = os.path.join(
        args.output_dir,
        f"average_top{args.num_taxa}_{args.treatment_column}_{args.filter_column}.tsv"
    )
    output_df.to_csv(filename, sep='\t', index=False)


# ---------------------------
# Per-treatment output
# ---------------------------
if args.mode == 'per_treatment':
    logging.info("Generating per-treatment files")

    per_treatment_dir = os.path.join(args.output_dir, "split_outputs")
    logging.info("Output files will be located in: %s", per_treatment_dir)

    os.makedirs(per_treatment_dir, exist_ok=True)

    for treatment in treatments:
        if treatment not in treatment_avg:
            continue

        safe_treatment = safe_name(treatment)

        top_taxa = treatment_avg[treatment].head(args.num_taxa).index.tolist()
        row = build_row(treatment_avg[treatment], top_taxa, treatment)
        df_out = pd.DataFrame([row])
        df_out = sort_columns_by_abundance(df_out)

        filename = os.path.join(
            per_treatment_dir,
            f"{args.treatment_column}_{safe_treatment}_{args.filter_column}.tsv"
        )

        df_out.to_csv(filename, sep='\t', index=False)


# ---------------------------
logging.info("Done.")