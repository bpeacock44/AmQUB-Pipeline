#!/usr/bin/env python3

import sys

# Define helper functions
def load_qiime_asv_table(fp):
    # Load a QIIME-happy text ASV table file
    with open(fp, 'r') as file:
        # Read all lines
        lines = file.readlines()

    # Extract column names from the first row
    col_names = lines[0].strip().split('\t')
    
    # Load data excluding the header
    data = [line.strip().split('\t') for line in lines[1:]]

    return col_names, data

def sort_qiime_asv_table(col_names, data, sortby="row", normalize_sort=False):
    # Sort tbl descending by rowSums or ascending by colSums
    numcols = [idx for idx, col_name in enumerate(col_names) if col_name.isdigit()]
    if sortby == "row":
        if normalize_sort:
            for row in data:
                row_sum = sum(float(row[col_idx]) for col_idx in numcols)
                if row_sum != 0:
                    for col_idx in numcols:
                        row[col_idx] = float(row[col_idx]) / row_sum
        data.sort(key=lambda x: sum(float(x[col_idx]) for col_idx in numcols), reverse=True)
    elif sortby == "col":
        data.sort(key=lambda x: [x[idx] for idx in numcols])
    else:
        raise ValueError("Values for 'sortby' must be either 'row' or 'col'")
    return col_names, data

def find_last_header_line(fp):
    # Your implementation of finding the last header line
    pass

def save_table(col_names, data, output_file):
    # Save the table to a new file
    with open(output_file, 'w') as file:
        # Write header
        file.write('\t'.join(col_names) + '\n')
        # Write data
        for row in data:
            file.write('\t'.join(row) + '\n')
    print(f"Sorted ASV table saved to {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python qiime_table_sorter.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Load ASV table
    col_names, data = load_qiime_asv_table(input_file)

    # Sort columns alphabetically
    col_names, data = sort_qiime_asv_table(col_names, data, sortby="col", normalize_sort=False)

    # Sort rows descending by rowSums
    col_names, data = sort_qiime_asv_table(col_names, data, sortby="row", normalize_sort=False)

    # Save the table to a new file
    save_table(col_names, data, output_file)

if __name__ == "__main__":
    main()
