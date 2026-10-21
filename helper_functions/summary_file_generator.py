#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse


def provided(path):
    return bool(path) and path.lower() != "null"


def process_data_and_write_summary(norm, raw, tax, tax_c, mixed, output_file):

    use_blast = provided(tax)
    use_classifier = provided(tax_c)
    use_mixed = provided(mixed)

    if not use_blast and not use_classifier:
        raise ValueError(
            "At least one taxonomy source must be provided."
        )

    # Read normalized abundance table
    df = pd.read_csv(norm, sep="\t", comment=None, header=0)
    df.rename(columns={df.columns[0]: "ID"}, inplace=True)

    # Remove taxonomy already embedded in the table
    df.drop(columns=["taxonomy"], inplace=True, errors="ignore")

    # Identify abundance columns before adding metadata
    metadata_cols = {"ID", "Sequence", "sequence"}
    abundance_cols = [
        col for col in df.columns
        if col not in metadata_cols
        and pd.api.types.is_numeric_dtype(df[col])
    ]

    # BLAST taxonomy
    if use_blast:
        tax_df = pd.read_csv(tax, sep="\t", header=0)
        tax_df.columns = [
            "ID", "taxonomy", "bitscore", "per_ID", "per_qcov"
        ]
        df = df.merge(
            tax_df[["ID", "taxonomy", "per_ID", "per_qcov"]],
            on="ID",
            how="left",
            validate="many_to_one"
        )
    else:
        df["taxonomy"] = np.nan
        df["per_ID"] = np.nan
        df["per_qcov"] = np.nan

    # Classifier taxonomy
    if use_classifier:
        classifier_df = pd.read_csv(tax_c, sep="\t", header=0)
        classifier_df.columns = [
            "ID", "c_taxonomy", "confidence"
        ]
        df = df.merge(
            classifier_df,
            on="ID",
            how="left",
            validate="many_to_one"
        )

    # Mean relative abundance
    df["avg_abun"] = df[abundance_cols].mean(
        axis=1, skipna=True
    )

    # Mixed-family analysis
    if use_mixed:
        top10_df = pd.read_csv(
            mixed, sep="\t", comment=None, header=None
        )
        df["mixed"] = np.where(
            df["ID"].isin(top10_df.iloc[:, 0]),
            "yes", "no"
        )

    # Explicit column ordering
    metadata_order = [
        "ID", "Sequence", "taxonomy", "per_ID", "per_qcov"
    ]

    if use_classifier:
        metadata_order += ["c_taxonomy", "confidence"]

    metadata_order.append("avg_abun")

    if use_mixed:
        metadata_order.append("mixed")

    # Preserve other columns, including abundance data
    remaining_cols = [
        col for col in df.columns
        if col not in metadata_order
    ]

    df = df[metadata_order + remaining_cols]

    # Read raw counts
    rawdf = pd.read_csv(raw, sep="\t", comment=None, header=0)

    rawdf.drop(
        columns=[
            rawdf.columns[0],
            "taxonomy",
            "Sequence",
            "sequence"
        ],
        inplace=True,
        errors="ignore"
    )

    # Calculate raw count totals
    sums = rawdf.sum(numeric_only=True)

    # Construct total-count row by column name
    total_row = {col: "nd" for col in df.columns}

    for col in ["per_ID", "per_qcov", "confidence", "avg_abun"]:
        if col in total_row:
            total_row[col] = 0

    for col, value in sums.items():
        if col in total_row:
            total_row[col] = value

    df = pd.concat(
        [pd.DataFrame([total_row]), df],
        ignore_index=True
    )

    # Save summary
    df.to_csv(
        output_file,
        sep="\t",
        index=False,
        quoting=False
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate detailed informational TU table."
    )

    parser.add_argument("norm")
    parser.add_argument("raw")
    parser.add_argument("tax")
    parser.add_argument("tax_c")
    parser.add_argument("mixed")
    parser.add_argument("output_file")

    args = parser.parse_args()

    process_data_and_write_summary(
        args.norm,
        args.raw,
        args.tax,
        args.tax_c,
        args.mixed,
        args.output_file
    )


if __name__ == "__main__":
    main()
