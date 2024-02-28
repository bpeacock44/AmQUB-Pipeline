#!/usr/bin/env python
'''
This is the 3rd script of the following 3 scripts:
1) check_barcode_collisions.pl
2) filter_barcode_noncollisions.py
3) fastq_convert_mm2pm_barcodes.py

which facilitate the recovery of reads that have barcodes not
perfectly matched to your expected/intended barcodes.

This script:
1) loads an output file of script 2)
2) loads the fastq barcode file corresponding to that output file
3) converts mismatched (MM) barcdoes in 2) into their
   perfect-match (PM) barcode versions

A WORD OF CAUTION:
When converting MM barcodes to PM barcodes, one must consider how "close"
(by Hamming distance) other true/perfect-match barcodes are to 
each other. If the distance between PM barcodes is small (eg., 3 or less) 
it is NOT ADVISABLE to convert all possible MM barcdes, as this could 
cause many reads from various samples to "bleed into" each other, and 
thus potentially confound your results.
'''

import argparse
import sys
import os
import re
from os.path import isfile
#import datetime


def read_mm2pm_bcfile(opts):
    '''Load the barcodes file'''
    
    #dictionary: d[mmbc] = pmbc
    bcs_d = {}
    
    #load barcodes file
    with open(opts.mm2pm_barcodes_fp, 'r') as mm2pm_bcfile:

        for line in mm2pm_bcfile:

            #Format example of the input file:
            # #MMBCs  PMBCs
            # TACACT  CACACT
            # CCCACT  CACACT
            # CGCACT  CACACT
            # CAGACT  CACACT
            # CACTCT  CACACT
            # ACTCTC  TCTCTC
            # TATCTC  TCTCTC
            # TGTCTC  TCTCTC
            # TCACTC  TCTCTC

            #capture values in the first three columns of (some) lines
            searchObj = re.search(r'^([ACTGN]+)\t([ACTGN]+)', line, re.M)

            #if it looks like we've got two barcodes:
            if searchObj and len(searchObj.groups()) == 2:
                mmbc = searchObj.group(1)
                pmbc = searchObj.group(2)
                
                #compare lengths as an additional check
                if len(mmbc) == len(pmbc):
                    #add to the dict
                    bcs_d[mmbc] = pmbc
                
    mm2pm_bcfile.close()

    return bcs_d
#

def read_convert_save(opts, bcs_d):
    '''read barcodes, convert, save result'''

    #notify
    print ("Saving  ["+opts.output_fp+"]")

    #open output file for writing
    fqfileout = open(opts.output_fp, 'w')

    changed_barcodes_count = 0
    total_barcodes_count = 0

    #load barcodes file
    with open(opts.input_fp, 'r') as fqfile:

        line_number = 0
        if opts.file_type == "barcode":
            for line in fqfile:

                #if we are on a barcode line...
                if line_number % 4 == 1:

                    total_barcodes_count += 1

                    #trim return char at the end
                    bc = line[:-1]

                    #convert it, if we should
                    if bc in bcs_d:
                        line = bcs_d[bc]+"\n"
                        changed_barcodes_count += 1

                #save it, converted or not
                fqfileout.write(line)
                line_number += 1

        elif opts.file_type == "read":
            for line in fqfile:
                
                #if we are on an @ line...
                if line_number % 4 == 0:

                    total_barcodes_count += 1

                    #trim return char at the end
                    elems = re.split(':', line.rstrip())
                    bc = elems[-1]
                    del elems[-1]

                    #convert it, if we should
                    if bc in bcs_d:
                        elems = elems + [bcs_d[bc]+"\n"]
                        line = ':'.join(elems)
                        changed_barcodes_count += 1

                #save it, converted or not
                fqfileout.write(line)
                line_number += 1

    print (str(total_barcodes_count)+" total barcodes in barcode file")
    print (str(changed_barcodes_count)+" MM barcodes were changed to PM barcodes")
    # pct_changed = changed_barcodes_count / total_barcodes_count * 100
    # print ("("+str("{0:.2f}".format(pct_changed))+"% of MM barcodes were changed to PM barcodes)")
    #TODO: count PM barcodes. pct_changed should calculate the % increase in PM barcodes, not all barcodes in the fastq file


def get_fastq_type(opts):
    '''get fastq type (currently detects either Illumina "RAW" or "USEARCH")'''

    #create a parser object
    parser = ParseFastQ(opts.input_fp)

    #read the first record
    for rec in parser:
        #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        #Example of RAW:
        #@M02457:244:000000000-BV7KW:1:1102:18589:1920 1:N:0:GTCCTGATGCCT
        #
        #Example of USEARCH:
        #@M02457:244:000000000-BV7KW:1:1102:15655:1826;sample=B023.81;
        #
        elems = re.split(':', rec[0])
        if len(elems) == 10:
            return "RAW"
        elif len(elems) == 7:
            return "USEARCH"
        else:
            print ("get_fastq_type: *** Error *** Fastq type not recognized (Illumina RAW or USEARCH")
            sys.exit(0)



def get_opts(argv):
    #create an ArgumentParser object
    parser = argparse.ArgumentParser(
        description="Convert selected mismatch barcodes in a fastq file into their corresponding "
                    +"perfect-match barcodes.")

    #add arguments
    requiredArg = parser.add_argument_group('required arguments')
    requiredArg.add_argument("-i", "--input_fp", help="The input fastq filepath with mismatch-barcodes "
                        +"that you want to change into perfect-match barcodes. It should be the same "
                        +"fastq file that [check_barcode_collisions.pl] was applied to.")
                        
    requiredArg.add_argument("-m", "--mm2pm_barcodes_fp", help="The filepath containing mismatch "
                        +"and perfect-match barcodes. A file of this type can be created by the "
                        +"script [filter_barcode_noncollisions.py]")

    requiredArg.add_argument("-o", "--output_fp", help="The output fastq filepath.")

    parser.add_argument("-t", "--file_type", help="Type of fastq file (barcode or read)")

    parser.add_argument("-f", "--force_overwrite", help="Overwrite the output fastq file if it exists.",
                        action='store_true')

    #parse the command line parameters
    opts = parser.parse_args()

    return opts


def check_opts(opts):
    #verify the given parameters are sane
    if not os.path.isfile(opts.input_fp):
        print ("** Error ** Input file ["+opts.input_fp+"] does not exist!")
        sys.exit(0)
    if not os.path.isfile(opts.mm2pm_barcodes_fp):
        print ("** Error ** Input file ["+opts.mm2pm_barcodes_fp+"] does not exist!")
        sys.exit(0)
    if opts.input_fp == opts.output_fp:
        print ("** Error ** Input and Output filenames must not be identical!")
        sys.exit(0)
    if os.path.isfile(opts.output_fp) and not opts.force_overwrite:
        print ("** Warning ** Output file ["+opts.output_fp+"] already exists!")
        sys.exit(0)
    if opts.file_type is None:
        print ("** Error ** Please supply a value for option [--file_type]!")
        print ("Choices are 'barcode' or 'read'")
        sys.exit(0)
    elif opts.file_type not in ['read','barcode']:
        print ("** Error ** Unknown file_type ["+str(opts.file_type)+"]!")
        print ("Choices are 'barcode' or 'read' only")
        sys.exit(0)



def main(argv):
    '''main function'''

    opts = get_opts(argv)

    check_opts(opts)

    print ("Loading ["+opts.mm2pm_barcodes_fp+"]")
    #load the barcodes into a dict
    bcs_d = read_mm2pm_bcfile(opts)

    print ("Loading ["+opts.input_fp+"]")
    #read the fastq file, convert barcodes, save to output file
    read_convert_save(opts, bcs_d)






if __name__ == "__main__":
    main(sys.argv)


