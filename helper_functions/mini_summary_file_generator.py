#!/usr/bin/env python3
"""Generate the Part 4 taxonomy summary; BLAST and classifier are independently optional."""
import argparse
import pandas as pd


def available(path):
    return bool(path) and path.lower() != 'null'


def read_fasta(path):
    records = []
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, f'>{header}#{"".join(chunks)}'))
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line)
        if header is not None:
            records.append((header, f'>{header}#{"".join(chunks)}'))
    return pd.DataFrame(records, columns=['ID', 'sequence'])


def process_data_and_write_summary(norm, fa, tax, tax_c, output_file):
    if not (available(tax) or available(tax_c)):
        raise ValueError('At least one taxonomy source (BLAST or classifier) is required.')

    df = pd.read_csv(norm, sep='\t', comment=None, header=1)
    df.rename(columns={df.columns[0]: 'ID'}, inplace=True)
    abundance_cols = df.select_dtypes(include='number').columns.tolist()

    if available(tax):
        blast = pd.read_csv(tax, sep='\t', comment=None, header=0)
        blast.columns = ['ID', 'taxonomy', 'blast_bitscore', 'blast_per_ID', 'blast_per_qcov']
        df = df.merge(blast[['ID', 'taxonomy', 'blast_per_ID', 'blast_per_qcov']], on='ID', how='left', validate='many_to_one')
    else:
        df['taxonomy'] = pd.NA
        df['blast_per_ID'] = pd.NA
        df['blast_per_qcov'] = pd.NA

    if available(tax_c):
        classifier = pd.read_csv(tax_c, sep='\t', comment=None, header=0)
        classifier.columns = ['ID', 'classifier_taxonomy', 'classifier_confidence']
        df = df.merge(classifier, on='ID', how='left', validate='many_to_one')

    df['avg_abun'] = df[abundance_cols].mean(axis=1, skipna=True)
    df.drop(columns=abundance_cols, inplace=True)
    df = df.merge(read_fasta(fa), on='ID', how='left', validate='many_to_one')
    df.to_csv(output_file, sep='\t', index=False, quoting=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('norm')
    parser.add_argument('fa')
    parser.add_argument('tax', help='BLAST taxonomy file or NULL')
    parser.add_argument('tax_c', help='Classifier taxonomy file or NULL')
    parser.add_argument('output_file')
    args = parser.parse_args()
    process_data_and_write_summary(args.norm, args.fa, args.tax, args.tax_c, args.output_file)


if __name__ == '__main__':
    main()
