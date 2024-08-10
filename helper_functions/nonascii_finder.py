#!/usr/bin/env python3

import os
import glob
import argparse

def has_non_ascii(line):
    for c in line:
        if not (ord(c) < 128 or c == '\t'):  # Check if character is ASCII or tab
            return True
    return False

def replace_non_ascii(input_file):
    temp_file = input_file + '.temp'
    replacements = []

    with open(input_file, 'r', encoding='utf-8') as fin, open(temp_file, 'w', encoding='utf-8') as fout:
        for line in fin:
            cleaned_line = ''.join(c if ord(c) < 128 else '?' for c in line)  # Replace non-ASCII with '?'
            if cleaned_line != line:
                replacements.append((line.strip(), cleaned_line.strip()))
            fout.write(cleaned_line + '\n')  # Write cleaned line to temp file

    os.replace(temp_file, input_file)
    return replacements

def main(pattern):
    files = glob.glob(pattern)
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return
    
    print(f"Processing files matching pattern: {pattern}")
    for filename in files:
        print(f"Processing file: {filename}")
        replacements = replace_non_ascii(filename)
        if replacements:
            for before, after in replacements:
                print(f"Replaced non-ASCII characters:\nBefore: {before}\nAfter : {after}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace non-ASCII characters in files matching a pattern.")
    parser.add_argument("pattern", type=str, help="Pattern to match files (e.g., '*final*/asvs/rep_set/assgntax/*seqs_chimera_filtered*')")
    args = parser.parse_args()

    main(args.pattern)
