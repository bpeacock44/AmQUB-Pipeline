#!/usr/bin/env python

'''
This script is for preparing an Illumina fastq read file for SRA upload,
removing reads who's barcodes or SampleIDs are NOT FOUND in the mapping file.

Currently this works for ONLY Illumina or UPARSE-FORMATED fastq headers

Example:
filter_fastx.py -i JB81_FC869A1P1.fastq -m JB81_map.txt -o sra/

Copyright 2019, Paul Ruegger
'''


from __future__ import print_function, division

import argparse
import sys
import os
import re
from os.path import isfile



class ParseFastQ(object):
    #Copied from: https://scipher.wordpress.com/2010/05/06/simple-python-fastq-parser/
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Usage:
        parser = ParseFastQ(filePath)

        parser.next()
        -OR-
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)


def read_sampleIDs(opts):
    if opts.fasta_fp:
        print ("Loading ["+opts.fasta_fp+"]")
        ids_d = read_fasta(opts)
    elif opts.map_fp:
        print ("Loading ["+opts.map_fp+"]")
        ids_d = read_map(opts)

    return ids_d


def read_fasta(opts):
    '''Read the fasta file'''
    #Format example:
    # >B001.81
    # CTCGACTACTGA

    with open(opts.fasta_fp, 'r') as idfile:
        lines = (line.rstrip() for line in idfile)
        lines = list(line for line in lines if line) # get only non-blank lines
    idfile.close()

    #get ids and barcodes
    ids_d = {}
    for i in xrange(0,len(lines),2):
        id = lines[i][1:]
        bc = lines[i+1]

        if bc in ids_d:
            print("*** WARNING *** Barcode file contains duplicate barcodes! Barcodes MUST be unique per map file!")
            print(bc)
            sys.exit(0)
        else:
            ids_d[bc] = id

    keys = ids_d.keys()
    print("  Loaded "+str(len(keys))+" sampleIDs. Eg: "+keys[1])
    return ids_d


def read_map(opts):
    '''Read the mapping file'''
    #Format example of the input file:
    # #SampleID	BarcodeSequence
    # B001.81	CTCGACTACTGA

    ids_d = {}
    with open(opts.map_fp, 'r') as idfile:
        for line in idfile:
            #skip header(s)
            if "#" in line[0:1]:
                continue

            #capture values in the first two columns
            searchObj = re.search(r'^(\S+)\t(\S+)', line, re.M)
            if searchObj and len(searchObj.groups()) == 2:
                id = searchObj.group(1)
                bc = searchObj.group(2)
                if bc in ids_d:
                    print("*** WARNING *** Mapping file contains duplicate barcodes! Barcodes MUST be unique per map file!")
                    sys.exit(0)
                else:
                    ids_d[bc] = id

    idfile.close()
    keys = ids_d.keys()
    print("  Loaded "+str(len(keys))+" sampleIDs. Eg: "+keys[0])
    return ids_d


def read_fastq_save_separate(opts, ids_d, fqtype):
    '''read fastq, convert, save result'''

    #notify
    print ("Saving  ["+opts.output_dir+"<id>.fastq]")

    #create output filenames and open one for each ID
    fqfilehandles = {}
    for bc in ids_d.keys():
        id = ids_d[bc]
        fp = opts.output_dir+id+".fastq"
        fh = open(fp, 'w')
        fqfilehandles[id] = fh

    #create a parser object
    parser = ParseFastQ(opts.input_fp)

    #we'll keep track of counts
    map_barcodes_count = 0
    total_barcodes_count = 0

    if fqtype == "RAW":
        #read the fastq file
        for rec in parser:
            total_barcodes_count += 1

            #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
            #@M02457:244:000000000-BV7KW:1:1102:15655:1826 1:N:0:ACCAGCCTATAG
            elems = re.split(':', rec[0])
            bc = elems[-1]

            #save this read if we should
            if bc in ids_d:
                id = ids_d[bc]
                map_barcodes_count += 1

                #write record to the correct file
                fqfilehandles[id].write('\n'.join(rec)+"\n")
                
    elif fqtype == "USEARCH":
        #read the fastq file
        for rec in parser:
            total_barcodes_count += 1

            #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
            #@M02457:244:000000000-BV7KW:1:1102:15655:1826;sample=B023.81;
            elems = re.split('sample=', rec[0])
            el = elems[-1]
            id = el[:-1]

            #save this read if we should
            if id in ids_d.values():
                map_barcodes_count += 1

                #write record to the correct file
                fqfilehandles[id].write('\n'.join(rec)+"\n")

    print (str(total_barcodes_count)+" total reads in fastq file")
    print (str(map_barcodes_count)+" reads with a map file SampleID were saved")


def read_fastq_save(opts, ids_d, fqtype):
    '''read fastq, convert, save result to a single file'''

    #notify
    print ("Saving  ["+opts.output_fp+"]")

    #create output file and open
    fqfileout = open(opts.output_fp, 'w')

    #create a parser object
    parser = ParseFastQ(opts.input_fp)

    #we'll keep track of counts
    map_barcodes_count = 0
    total_barcodes_count = 0

    if fqtype == "RAW":
        #read the fastq file
        for rec in parser:
            total_barcodes_count += 1

            #rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
            #@M02457:244:000000000-BV7KW:1:1102:15655:1826 1:N:0:ACCAGCCTATAG
            elems = re.split(':', rec[0])
            bc = elems[-1]

            #save this read if we should
            if bc in ids_d:
                map_barcodes_count += 1

                #write record to the correct file
                fqfileout.write('\n'.join(rec)+"\n")
        
    else:
        print ("*** Error *** Fastq type not recognized (must be Illumina RAW)")
        sys.exit(0)


    print (str(total_barcodes_count)+" total reads in fastq file")
    print (str(map_barcodes_count)+" reads with a map file SampleID were saved")


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
        description="Separate a fastq file into separate files by sampleIDs in a mapping file.")

    #add arguments
    requiredArg = parser.add_argument_group('required arguments')
    requiredArg.add_argument("-i", "--input_fp", help="The input fastq filepath.")

    parser.add_argument("-m", "--map_fp", help="A mapping file containing #SampleID<tab>BarcodeSequence columns.")
    parser.add_argument("-b", "--fasta_fp", help="A fasta file containing >SAMPLEID<cr>BARCODE rows.")

    parser.add_argument("-d", "--output_dir", help="The output directory to save fastq files.")
    parser.add_argument("-o", "--output_fp", help="The output fastq file.")

    parser.add_argument("-f", "--force_overwrite", help="Overwrite the output fastq file if it exists.",
                        action='store_true')

    #parse the command line parameters
    opts = parser.parse_args()

    return opts


def check_opts(opts):
    if not opts.input_fp:
        print ("** Error ** Please specify a value for --input_fp!")
        sys.exit(0)
    if not opts.map_fp and not opts.fasta_fp:
        print ("** Error ** Please specify a value for --map_fp OR --fasta_fp!")
        sys.exit(0)
    if opts.map_fp and opts.fasta_fp:
        print ("** Error ** Cannot use both --map_fp AND --fasta_fp options!")
        sys.exit(0)
    #set a single var for the barcode filepath
    if opts.map_fp:
        barcode_fp = opts.map_fp
    if opts.fasta_fp:
        barcode_fp = opts.fasta_fp


    #verify the given parameters are sane
    if not os.path.isfile(opts.input_fp):
        print ("** Error ** Input file ["+opts.input_fp+"] does not exist!")
        sys.exit(0)
    if not os.path.isfile(barcode_fp):
        print ("** Error ** Input file ["+barcode_fp+"] does not exist!")
        sys.exit(0)
    if opts.input_fp == opts.output_dir:
        print ("** Error ** Input and Output filenames must not be identical!")
        sys.exit(0)

    if not opts.output_dir and not opts.output_fp:
        print ("** Error ** Please specify either --output_dir OR --output_fp!")
        sys.exit(0)
    if opts.output_dir and opts.output_fp:
        print ("** Error ** Please specify either --output_dir OR --output_fp!")
        sys.exit(0)


    if opts.output_dir:
        #add '/' if needed
        if(opts.output_dir[-1] != "/"):
            opts.output_dir += "/"

        if os.path.isdir(opts.output_dir) and not opts.force_overwrite:
            print ("** Warning ** Output directory ["+opts.output_dir+"] already exists!")
            sys.exit(0)
        elif not os.path.exists(opts.output_dir):
            os.makedirs(opts.output_dir)
            print ("Created output dir ["+opts.output_dir+"]")
    
    if opts.output_fp:
        if os.path.isfile(opts.output_fp) and not opts.force_overwrite:
            print ("** Warning ** Output file ["+opts.output_fp+"] already exists!")
            sys.exit(0)
        elif not os.path.exists(opts.output_fp):
            print ("Saving to ["+opts.output_fp+"]")

    # if opts.file_type is None:
        # opts.file_type = "barcode"
    # elif opts.file_type not in ['read','barcode']:
        # print ("** Error ** Unknown file_type ["+str(opts.file_type)+"]!")
        # print ("Choices are 'barcode' or 'read' only")
        # sys.exit(0)


def main(argv):
    '''main function'''

    opts = get_opts(argv)

    check_opts(opts)

    #determine the input fastq's format (Illumina 'RAW' or 'USEARCH')
    fqtype = get_fastq_type(opts)

    #load the map's SampleIDs and BarcodeSequence's into a dict
    ids_d = read_sampleIDs(opts)

    #read the input fastq file, separate by BarcodeSequence's or SampleIDs, save to output file(s)
    if opts.output_dir:
        print ("Loading ["+opts.input_fp+"] and saving output to ["+opts.output_dir+"]")
        read_fastq_save_separate(opts, ids_d, fqtype)
    elif opts.output_fp:
        read_fastq_save(opts, ids_d, fqtype)
    else:
        print ("Something is rotten in Denmark")




if __name__ == "__main__":
    main(sys.argv)