#!/usr/bin/env python
'''
This script creates a list of mismatched fastq barcodes and their corresponding perfect-match 
barcodes from the output file of [check_barcode_collisions.pl]. The output format is:

<mismatched barcode><tab><barcode>
<mismatched barcode><tab><barcode>
...

By default the script will discard mismatched barcodes if their counts in the fastq file are 
greater than the count of their corresponding perfect-match barcode, as this may indicate that
mismatches from multiple barcodes have contributed to its count. This behavior can be overridden
by using the -k switch.

The output of this script can be used as input to the [fastq_convert_mm2pm_barcodes.py] script, 
which creates a new fastq barcode file whereby all acceptable mismatched barcodes have been 
changed to their perfect match versions.

This script is step 2 of 3 scripts:
1) check_barcode_collisions.pl
2) filter_barcode_noncollisions.py
3) fastq_convert_mm2pm_barcodes.py

'''

#pylint $(which filter_barcode_noncollisions.py) | less
#filter_barcode_noncollisions.py -i uFQBC_PR04_FC233L1P3.fastq_BC6_M2.txt | less

from __future__ import print_function, division

import argparse
import sys
import os
import re
from os.path import isfile
import datetime
import csv


class Barcode:
    '''The Barcode class'''
    def __init__(self, seq, fqcount):
        self.seq = seq
        self.fastq_count = int(fqcount)
        self.fastq_count_with_mm_barcodes = self.fastq_count
        self.all_fastq_counts = self.fastq_count

        self.mm_count = None
        self.assigned_sampleID = None
        self.is_pm = False
        self.mmbcs = [] #holds barcode objects with mismatches to itself

    def get_seq(self):
        return self.seq

    def get_fastq_count(self):
        return self.fastq_count

    def get_fastq_count_with_mm_barcodes(self):
        return self.fastq_count_with_mm_barcodes

    def set_mm_count(self, mms):
        self.mm_count = mms
        if self.mm_count == 0:
            self.is_pm = True

    def get_mm_count(self):
        return int(self.mm_count)


    def set_assigned_SampleID(self, sampleID):
        self.assigned_sampleID = str(sampleID)

    def get_assigned_SampleID(self):
        assert (self.assigned_sampleID is not None), "'assigned_sampleID' has not been set!"
        return self.assigned_sampleID

    def add_mmbc(self, mmbc):
        assert (self.is_pm is True), "Adding barcodes to barcodes with mismatches is not allowed!"
        self.mmbcs.append(mmbc)
        self.fastq_count_with_mm_barcodes += mmbc.get_fastq_count()

    def add_to_all_fastq_counts(self, count):
        self.all_fastq_counts += int(count)
        
    def get_all_fastq_counts(self):
        return self.all_fastq_counts

    def get_mmbcs(self):
        #return sorted so diff checks in different runs show real diffs (not diffs due to order)
        return sorted(self.mmbcs)
#


def read_uFQBC_file(opts):
    '''
    Reads the input file and returns a list of Perfect-Match barcodes, where each
    PM barcode may contain a list of mismatch barcodes to itself, and also a
    total count of the MM barcodes.
    '''
    pm_Barcodes = []

    #load file
    with open(opts.input_fp, 'r') as infile:

        current_PM_barcode = None
        for line in infile:

            #Format example of the input file (output of [check_barcode_collisions.pl]):
            # bc	CACACT	239	[F01]
            # m1	AACACT	448	0
            # *1	GACACT	181	0	[ F88:GAAACT F01:CACACT ]
            # m1	TACACT	90	0
            # *1	CTCACT	81	1	[ F27:CTCAAT F01:CACACT ]
            # m1	CCCACT	64	1

            #capture values in the first three columns of (some) lines
            searchObj = re.search(r'^(bc|m['+str(opts.mismatches)+'])\t(\w+)\t(\d+)', line, re.M)

            #keep matched info
            if searchObj and len(searchObj.groups()) == 3:

                #clarify meaning of captured groups
                bctype  = searchObj.group(1) #the type of barcode (eg: bc=perfect match, m1=mismatch)
                bcseq   = searchObj.group(2) #the barcode or mismatched barcode sequence
                fqcount = searchObj.group(3) #the count of this barcode in the fastq file

                #initialize a Barcode object
                bc = Barcode(bcseq, fqcount)

                #if bc is a perfect-match barcode
                if bctype == 'bc':
                    #change its name so we won't lose the reference to it
                    current_PM_barcode = bc
                    #set its mismatch count to zero
                    current_PM_barcode.set_mm_count(0)

                    #set the SampleID it represents
                    searchObj2 = re.search(r'\[(\S+)\]', line, re.M)
                    if searchObj2:
                        current_PM_barcode.set_assigned_SampleID(searchObj2.group(1))

                    #add to the list of perfect match barcodes
                    pm_Barcodes.append(current_PM_barcode)

                #if bc is a mismatched barcode
                elif bctype.startswith('m'):
                    #set mismatch count
                    bc.set_mm_count(bctype[1])
                    
                    #TODO: Store ALL mm barcodes in Barcode but include/set a flag to denote
                    # whether or not they have passed the (arbitrary) filtering criteria?
                    
                    if opts.keep_larger_mismatches:
                        current_PM_barcode.add_mmbc(bc)
                    #add bc to its perfect match Barcode list IF it meets criteria
                    elif bc.get_fastq_count() < current_PM_barcode.get_fastq_count():
                        current_PM_barcode.add_mmbc(bc)

                    #this is a hack to keep track of counts without the filtering above
                    current_PM_barcode.add_to_all_fastq_counts(fqcount)

    infile.close()

    return pm_Barcodes


def print_for_fastq_convert(opts, pm_Barcodes):
    '''
    Prints output in this format:
    
    #MISMATCHES     PERFECTMATCHES(ie., the actual barcodes)
    <mmbarcode1><tab><pmbarcode1>
    <mmbarcode2><tab><pmbarcode1>
    <mmbarcode3><tab><pmbarcode1>
    ...
    <mmbarcode7><tab><pmbarcode2>
    <mmbarcode8><tab><pmbarcode2>
    <mmbarcode9><tab><pmbarcode2>
    ...
    '''

    #get percent that mm counts contribute#
    
    all_pm_counts = 0
    pm_and_mm_counts = 0
    for pmbc in pm_Barcodes:
        all_pm_counts += pmbc.get_fastq_count()
        pm_and_mm_counts += pmbc.get_fastq_count_with_mm_barcodes()
    
    #calc percent
    pct = "{0:.1f}%".format((pm_and_mm_counts / all_pm_counts - 1) * 100)
    

    #print header info
    print ("# Command: ", os.path.basename(sys.argv[0]), ' '.join(sys.argv[1:]))
    print ("# Date:    ", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ("# PM Count:", all_pm_counts)
    print ("# +MMs:    ", pm_and_mm_counts)
    print ("# +MM%:    ", pct)
    
    print ("# Use this file as input to [fastq_convert_mm2pm_barcodes.py]")
    print ("# MMBCs=MisMatch BarCodes, PMBCs=PerfectMatch BarCodes")
    print ("#MMBCs\tPMBCs")

    for pmbc in pm_Barcodes:
        #get a list of this barcode's family of mismatched barcodes
        mmbcs = pmbc.get_mmbcs()
        
        # if pmbc.get_seq() == "CTCGACTACTGA":
            # print ("pm = CTCGACTACTGA, type(mmbcs) =", type(mmbcs), "len =", len(mmbcs))
            # # sys.exit(0)
        # if mmbcs is None:
            # print (pmbc.get_seq(), "may not have any mm barcodes")
            
        #print each one, and the pm to which it belongs
        for mmbc in mmbcs:
            print ('\t'.join([mmbc.get_seq(), pmbc.get_seq()]))




def print_for_viewing(opts, pm_Barcodes):

    # #hack to print barcodes from just these samples
    # pr04 = ['F01','F02','F03','F08','F09','F10','F11','F12','F13','F14','F15','F16','F17','F18','F19','F20','F21','F22','F23','F24','F25','F26','F27','F30','F31','F32','F33','F34','F35','F36','F37','F38','F39','F40','F41','F42','F43','F44','F45','F46','F47','F48','F49','F50','F51','F52','F53','F54','F55','F56','F58','F59','F60','F61','F62','F63','F64','F65','F66','F67','F68','F69','F70','F71','F72','F73']


    #print header info
    print ("# Command:", os.path.basename(sys.argv[0]), ' '.join(sys.argv[1:]))
    print ("# Date:   ", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ()

    #print results
    for pmbc in pm_Barcodes:

        #perfect match Barcode
        print ('\t'.join(["bc", pmbc.get_seq(), str(pmbc.get_fastq_count()), '['+pmbc.get_assigned_SampleID()+']']))

        #TODO?: add an option to define a list/subset of SampleIDs for printing (vs printing all)
        # if pmbc.get_assigned_SampleID() in pr04:
            # print ('\t'.join(["bc", pmbc.get_seq(), str(pmbc.get_fastq_count()), '['+pmbc.get_assigned_SampleID()+']']))
        # else:
            # continue


        #get new total read counts if mismatch barcodes are converted to perfect match barcodes
        newtot = pmbc.get_fastq_count_with_mm_barcodes()

        #a list of this barcode's family of mismatched barcodes
        mmbcs = pmbc.get_mmbcs()

        #format the mmbcs info
        mm_info = []
        for mmbc in mmbcs:
            mm_info.append(['m'+str(mmbc.get_mm_count()), mmbc.get_seq(), str(mmbc.get_fastq_count())])

        #add new total read counts to the last mmbc
        if len(mm_info) > 0:
            # mm_info[-1].append("(new total: "+str(newtot)+")")
            mm_info[-1].append(str(pmbc.get_fastq_count())+"->"+str(newtot)+' ('+str(pmbc.get_all_fastq_counts())+' if no filtering)')


        #print mm_info for each mmbc
        for mmbc in mm_info:
            print ('\t'.join(mmbc))

        print()



#TODO: implement this option (probably should also implement a default output filename)
def save_to_csv(dict, output_fp):
    w = csv.writer(open(output_fp, "w"))
    for key, val in dict.items():
        w.writerow([key, val])
# #read from csv:
# dict = {}
# for key, val in csv.reader(open("input.csv")):
    # dict[key] = val



def main(argv):
    '''main function'''
    #create an ArgumentParser object
    parser = argparse.ArgumentParser(
        description="Find non-colliding barcodes from a [check_barcode_collisions.pl] output "
        +"file. By default, this script will discard any mismatched barcode if its count in the "
        +"fastq file is greater than the count of its corresponding perfect match barcode.")

    #add arguments
    
    requiredArg = parser.add_argument_group('required arguments')
    requiredArg.add_argument("-i", "--input_fp", help="The input fastq filepath", required=True)
    requiredArg.add_argument("-m", "--mismatches", help="The mismatches to keep. Eg. '1' for 1mm,"
                        +" or '12' to consider both 1 and 2 bp mismatches", type=int, required=True)

    parser.add_argument("-k", "--keep_larger_mismatches", help="keep mismatch barcodes even when "+
                        "their counts are greater than their corresponding perfect-match barcode "+
                        "(Default: False)", action="store_true", required=False)
    
    MEgroup = parser.add_mutually_exclusive_group(required=True)
    MEgroup.add_argument("--output_for_viewing", help="The output will be formatted for viewing only."
                        +" Cannot be used with OUTPUT_FOR_FASTQ_CONVERT", action="store_true")
    MEgroup.add_argument("--output_for_fastq_convert", help="The output will be formatted for use "
                        +"as input into the [fastq_convert_mm2pm_barcodes.py.py] script.\n"
                        +"Cannot be used with OUTPUT_FOR_VIEWING", action="store_true")

    #parse the command line parameters
    opts = parser.parse_args()


    #verify that input file exists
    if not os.path.isfile(opts.input_fp):
        print ("** Error ** Input file ["+opts.input_fp+"] does not exist!")
        sys.exit(0)

    # #check for exclusive options (argparser checks this via 'add_mutually_exclusive_group')
    # if opts.output_for_viewing and opts.output_for_fastq_convert:
        # print ("** Option Error ** Options [OUTPUT_FOR_VIEWING] and [OUTPUT_FOR_FASTQ_CONVERT]"
               # +"cannot be used simultaneously")

    #notifiy if output_for_viewing
    if opts.output_for_viewing:
        print ("Printing for viewing", file=sys.stderr)
    #notifiy if output_for_fastq_convert (using fastq_convert_mm2pm_barcodes.py)
    elif opts.output_for_fastq_convert:
        print ("Printing for input into [fastq_convert_mm2pm_barcodes.py]", file=sys.stderr)


    #list to hold perfect match Barcodes (each with their internal list of mismatched Barcodes)
    pm_Barcodes = read_uFQBC_file(opts)


    #print results for viewing
    if opts.output_for_viewing:
        # print ("Printing for viewing", file=sys.stderr)
        print_for_viewing(opts, pm_Barcodes)

    #print results useful for converting fastq barcodes with fastq_convert_mm2pm_barcodes.py
    if opts.output_for_fastq_convert:
        # print ("Printing for input into [fastq_convert_mm2pm_barcodes.py]", file=sys.stderr)
        print_for_fastq_convert(opts, pm_Barcodes)



if __name__ == "__main__":
    main(sys.argv)