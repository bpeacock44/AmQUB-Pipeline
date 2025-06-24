#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
import datetime

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

def main():
    parser = argparse.ArgumentParser(description="Process ASV table and mapping file to generate average abundance data.")
    parser.add_argument('--file_path', required=True)
    parser.add_argument('--map_file', required=True)
    parser.add_argument('--column', required=True)
    parser.add_argument('--output_avg_abundance_directory', required=True)
    parser.add_argument('--output_prism_file', required=True)
    parser.add_argument('--num_taxa', type=int, default=3)
    args = parser.parse_args()

    setup_logger(args.output_avg_abundance_directory)

    try:
        logging.info(f"Reading ASV table from: {args.file_path}")
        with open(args.file_path) as f:
            first_line = f.readline().strip()
        skip_rows = 1 if first_line == "# Constructed from biom file" else 0

        df = pd.read_csv(args.file_path, sep='\t', skiprows=skip_rows, comment='~')
        df.rename(columns={df.columns[0]: "ASV_ID"}, inplace=True)

        if 'taxonomy' in df.columns:
            df['taxonomy'] = df['taxonomy'] + '_' + df['ASV_ID']
            df.drop(columns=['ASV_ID'], inplace=True)
        else:
            df.rename(columns={'ASV_ID': 'taxonomy'}, inplace=True)

        logging.info(f"Reading mapping file from: {args.map_file}")
        map_df = pd.read_csv(args.map_file, sep='\t', comment='~')
        map_df.rename(columns={'#SampleID': 'sampleID'}, inplace=True)

        if args.column not in map_df.columns:
            logging.error(f"Column '{args.column}' not found in mapping file.")
            exit(1)

        map_df = map_df[[args.column, 'sampleID']]
        treatments = map_df[args.column].dropna().unique()

        os.makedirs(args.output_avg_abundance_directory, exist_ok=True)

        df_transposed = df.set_index('taxonomy').T.reset_index().rename(columns={'index': 'sampleID'})
        merged_df = pd.merge(map_df, df_transposed, on='sampleID')
        merged_df = merged_df.dropna(subset=[args.column])

        treatment_top_taxa_dict = {}
        unique_top_taxa = set()
        treatment_avg_abundance_dict = {}

        for treatment in treatments:
            logging.info(f"Processing treatment: {treatment}")
            t_df = merged_df[merged_df[args.column] == treatment]
            abund = t_df.drop(columns=['sampleID', args.column])
            avg_abund = abund.mean().reset_index()
            avg_abund.columns = ['taxonomy', 'average_abundance']
            top_taxa = avg_abund.sort_values(by='average_abundance', ascending=False).head(args.num_taxa)['taxonomy'].tolist()
            treatment_top_taxa_dict[treatment] = top_taxa
            unique_top_taxa.update(top_taxa)

        final_top_taxa = sorted(unique_top_taxa)

        for treatment in treatments:
            t_df = merged_df[merged_df[args.column] == treatment]
            abund = t_df.drop(columns=['sampleID', args.column])
            avg_abund = abund.mean().reset_index()
            avg_abund.columns = ['taxonomy', 'average_abundance']
            treatment_avg_abundance_dict[treatment] = avg_abund.set_index('taxonomy')

        columns = ['treatment'] + final_top_taxa + ['other']
        rows = []

        for treatment in treatments:
            avg_abund = treatment_avg_abundance_dict[treatment]['average_abundance']
            row = {'treatment': treatment}
            row.update({otu: avg_abund.get(otu, 0) for otu in final_top_taxa})
            row['other'] = avg_abund.drop(index=final_top_taxa, errors='ignore').sum()
            rows.append(row)

        output_df = pd.DataFrame(rows, columns=columns)

        treatment_means = pd.DataFrame([
            treatment_avg_abundance_dict[treatment]['average_abundance']
            for treatment in treatments
        ]).fillna(0)

        overall_row = {'treatment': 'avg_abun'}
        overall_row.update({otu: treatment_means[otu].mean() for otu in final_top_taxa})
        overall_row['other'] = treatment_means.drop(columns=final_top_taxa, errors='ignore').sum(axis=1).mean()

        output_df = pd.concat([pd.DataFrame([overall_row]), output_df], ignore_index=True)

        numeric_cols = output_df.columns.difference(['treatment'])
        output_df[numeric_cols] = output_df[numeric_cols].div(output_df[numeric_cols].sum(axis=1), axis=0) * 100

        logging.info(f"Saving the output PRISM file to: {args.output_prism_file}")
        output_df.to_csv(args.output_prism_file, sep='\t', index=False)
        logging.info("Processing complete. Results saved.")

    except Exception as e:
        logging.exception("An error occurred during processing.")
        exit(1)

if __name__ == '__main__':
    main()
