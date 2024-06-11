#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print("Usage: ./qiime_table_sorter.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the data from the input file and store it in a list of lists
data = []
with open(input_file, 'r') as file:
    for line in file:
        data.append(line.strip().split('\t'))

# Extract headers and remove them from the data
headers = data[0]
data = data[1:]

# Convert numerical values from strings to integers (assuming all columns after the first are numeric)
for row in data:
    for i in range(1, len(row)):
        row[i] = int(row[i])

# Calculate row sums and add them as an additional column
for row in data:
    row.append(sum(row[1:]))

# Sort the data by the last column (row sums) in descending order
sorted_data = sorted(data, key=lambda x: x[-1], reverse=True)

# Join sorted data with headers
sorted_with_headers = [headers] + sorted_data

# Transpose the data to make it easier to work with columns
transposed_data = list(zip(*sorted_with_headers))

# Sort columns alphabetically (excluding the first column which contains IDs)
sorted_columns = sorted(transposed_data[1:], key=lambda x: x[0])

# Transpose back the sorted columns
sorted_with_sorted_columns = [transposed_data[0]] + list(zip(*sorted_columns))

# Store sorted_with_sorted_columns[0] as a separate list called asvs
asvs = sorted_with_sorted_columns[0]

# Remove sorted_with_sorted_columns[0]
del sorted_with_sorted_columns[0]

sorted_with_sorted_columns = [list(row) for row in sorted_with_sorted_columns]

for i in range(len(sorted_with_sorted_columns)):
    sorted_with_sorted_columns[i].insert(0, asvs[i])

# Write sorted_with_sorted_columns to the output file
with open(output_file, 'w') as file:
    for row in sorted_with_sorted_columns:
        file.write('\t'.join(map(str, row)) + '\n')
