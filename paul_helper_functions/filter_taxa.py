#!/usr/bin/env python3

import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Filter lines from a file based on specified strings.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input file")
    parser.add_argument('-o', '--output', required=True, help="Path to the output file")
    parser.add_argument('-r', '--remove', required=True, help="Comma-delimited strings to filter out. Example: 'string1,string2,string3'")
    #parser.add_argument('-h', '--help', action='store_true', help="Show this help message and exit")
    
    # Parse arguments
    args = parser.parse_args()

    # Split the remove argument into a list of strings
    remove_strings = args.remove.split(',')

    # Read input file and filter lines
    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            if not any(remove_str in line for remove_str in remove_strings):
                outfile.write(line)

if __name__ == "__main__":
    main()
