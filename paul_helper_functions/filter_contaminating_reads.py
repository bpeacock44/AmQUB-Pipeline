#!/usr/bin/env python3

import argparse
import sys
import os
import re
from collections import defaultdict
import json #for converting defaultdict to dict
from os.path import isfile
import pprint


#phiX taxonomic ID
phiX_txid = '10847'
singletRaccns = {}
singletKaccns = {}


def usage():
    x="""
Create lists of non-contaminating (target) and contaminating (non-target) OTUs by parsing local BLAST results of your sequences against the nt database

Usage: perl $0 [OPTIONS]

   [REQUIRED]
   -i <BLAST output>   # BLAST output (using the format options given in #3 below)
   -k <taxID file>     # TaxonIDs of intended target organisms that should be KEPT
   -e <taxID file>     # TaxonIDs of Environmental samples of intended target organisms (provisionally KEPT. see below)
   -r <taxID file>     # TaxonIDs of likely contaminant organisms that should be REJECTED
   -t <directory>      # Directory that contains or will contain AccnsWithDubiousTaxAssigns.txt. Typically the location USEARCH is located.
   
   [OPTIONAL]
   -m <merge.dmp>      # Load the merge.dmp file from NCBI to update taxIDs **Important if your local BLAST db is old but your taxIDs are new**
                          merge.dmp can be downloaded from: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
   -p <prefix>         # Prefix for output files. Useful for version designations if trying multiple filtering variations (eg., 'v1', etc.)
   -n <integer>        # Stop processing BLAST output after n hits. Useful for debugging
   -x <regex>          # Report results by groups found in OTU names with Regular Expression <regex>
                          Special use cases only. Eg.,  Velvet or PFOR2 contigs that have modified fasta names containing SampleIDs
                          The regex should find the SampleID and contig lengths. Eg: 'denovoA10contig.14_0_149'=~/denovo(\w\d+)\S+_(\d+)/
   -N                  # Reject all OTUs with [N]o BLAST hits. See related option -d. (Default: Keep)
   -f                  # Force overwrite of output files (including log)
   -v                  # Print version and exit
   -h                  # This help

IMPORTANT: Before using this script you must download/create information from NCBI as follows:
1) Search NCBI's taxonomy db for the taxonomic group you're targeting in your assay, excluding environmental samples
   eg.: txid2[subtree] NOT \"environmental samples\"[subtree]
2) Search NCBI for the environmental samples of the taxonomic group you're targeting in your assay
   eg.: txid2[subtree] AND \"environmental samples\"[subtree]
3) Search NCBI for the taxonomic group of any likely contaminants, excluding environmental samples (eg., host plant)
   eg.: Mus[subtree] NOT \"environmental samples\"[subtree]
** For each of steps 1-3, download the associated TaxonID numbers and save to file (eg., Bacteria__k_txid2_NOT_Environmental_Samples.txt)
4) Perform a local blastn of your sequences against the nt database with the following format parameters:
    OPTS=\"qseqid sseqid pident length mismatch evalue bit-score staxids stitle\"
    blastn -db nt -query repseqs.fna -out repseqs.fna.blastout -max_target_seqs 150 -evalue 0.001 -outfmt 7 \"\$OPTS\"
    NOTE: max_target_seqs can be adjusted but a minimum of 100 is recommended to ensure most or all top hits in nt are returned
    IMPORTANT: In order to find taxonIDs (staxids), blastn MUST have access to the files 'taxdb.btd' and 'taxdb.bti'. These files comprise
               NCBI's taxonID database, and they must be located in the same directory as your local NCBI nt database files
5) Run this script. Example command:
    filter_contaminating_reads_by_NCBI_blasts.pl -i repseqs.fna.blastout \\
    -k Bacteria__k_txid2_NOT_Environmental_Samples.txt,Archaea__k_txid2157_NOT_Environmental_Samples.txt \\
    -e Bacteria__k_txid2_AND_Environmental_Samples.txt,Archaea__k_txid2157_AND_Environmental_Samples.txt \\
    -r Mouse__f_txid10066_NOT_Environmental_Samples.txt
    -t ~/path/to/directory
   The script will output 4 files: [otus2keep.txt], [otus2reject.txt], [otus2withNohits.txt] and [otus2filter.log]
6) Use QIIME filtering scripts to remove/keep the contaminating/non-contaminating reads from your files. Example command:
    filter_otus_from_otu_table.py -i original_otu_table.biom -o filtered_otu_table.biom --otu_ids_to_exclude_fp otus2reject.txt

** IMPORTANT: Providing a new <merge.dmp> file SHOULD BE CONSIDERED MANDATORY if the taxIDs files are newer than the BLAST db **
   LIKEWISE, ensure that the downloaded taxIDs ARE NOT OLDER than the BLAST db
   NCBI's taxonomy is not curated. To prevent sequences with dubiously assigned taxonomy from confounding your filtering results, add
    their ACCN.ver IDs to a file named [AccnsWithDubiousTaxAssigns.txt] in the same directory as the TaxonID files
   Lastly, note that OTUs with blast hits to the PhiX taxonID (10847) are automatically rejected
"""
    print(x)
#


class ParseBlastout(object):
    '''Returns a blasthit-by-blasthit BLASTN parser analogous to file.readline()'''
    '''Modified from: https://scipher.wordpress.com/2010/05/06/simple-python-fastq-parser'''
    def __init__(self,filePath,headerSymbol="# BLASTN"):
        '''
        -Usage-
        parser = ParseBlastout(filePath)
        parser.next()
        -OR-
        for rec in parser:
            ... do something with rec ...
        #(rec is a list)
        '''
        self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSym = headerSymbol

    def __iter__(self):
        return self

    def __next__(self):
        '''Reads in next element. Returns: list'''
        blast_hit = []

        #read first line at current position
        last_pos = self._file.tell()
        line = self._file.readline()
        
        #ensure we're at a new blast_hit or we stop (probably at eof)
        if self._hdSym in line:
            self._currentLineNumber += 1
        elif len(line) < 2:
            raise StopIteration
        elif len(line) > 0:
            print(line)
            print("Something is wrong with this blastout file")
            sys.exit(0)
        else:
            raise StopIteration
            
        #read lines until we get to the start of a new blast_hit
        while True:
            last_pos = self._file.tell()
            line = self._file.readline()
            
            #are we at the start of a new blast_hit?
            if self._hdSym in line:
                #set read pointer back one line (so we'll read this line again next time)
                self._file.seek(last_pos)
                break
            #are we at the end of the file?
            elif "BLAST processed" in line:
                break
            elif line:
                self._currentLineNumber += 1
                blast_hit.append(line.strip('\n'))
            else:
                #we are at the end of the file
                break
                
        return blast_hit
#


class blastHit:
    '''A blastHit object contains all relevant information from
       a blast hit in a condensed, sorted, and accessible way'''
    def __init__(self, otu, size, hitlines, usersTaxIDs):
        #bitscore, pident, taxid
        #determine: count of bitscore/taxid pairs
        self.otu = otu
        self.size = size
        self.num_hits = len(hitlines)
        
        self.choice = None ### holds classification decision choice # temp/debugging
        
        #dictionaries (d[staxid][bitscore]...)
        self.E = {}
        self.K = {}
        self.R = {}
        self.U = {}
        self.N = {}
        self.init_des = []
        self.final_des = []
        
        #lists of taxids
        self.txidsE = []
        self.txidsK = []
        self.txidsR = []
        self.txidsU = []
        self.txidsN = []
        
        #temp vars
        E = []
        K = []
        R = []
        U = []
        init_des = {}
        
        #parse hitlines
        for hit in hitlines:
            #hit example:
            # ['LR031876.1', '97.368', '4.61e-27', '130', '3712', 'Brassica oleracea HDEM genome, scaffold: C7', 100]
            #which are: [accnver, pident, evalue, bitscore, staxid, stitle, qcovs]
            
            #get 1st staxid if more than one
            if ";" in hit[4]:
                staxid = hit[4].split(';')[0]
            else:
                staxid = hit[4]
                
            try:
                staxid = int(staxid)
            except:
                #skip this hit - it's an 'N/A' or some other invalid taxid
                continue
                
            #replace original
            hit[4] = staxid
            
            #classify taxids in this hit by comparing them to the user-supplied E/K/R taxid files
            if staxid not in usersTaxIDs:
                #put it in the "Unknown" pile
                U.append(hit)
                init_des["U"] = 0
            else:
                if usersTaxIDs[staxid] == "E":
                    E.append(hit)
                    init_des["E"] = 0
                elif usersTaxIDs[staxid] == "K":
                    K.append(hit)
                    init_des["K"] = 0
                elif usersTaxIDs[staxid] == "R":
                    R.append(hit)
                    init_des["R"] = 0
                else:
                    print("Could not find an E, K, R or U designator for staxid:", staxid)
                    sys.exit(0)
                    
        #handle case where there are no hits
        if self.num_hits == 0:
            init_des["N"] = 0
            self.N = {'NoTaxID': 0}
            
        #sort all designations for consistent spelling
        self.init_des = sorted(init_des.keys())
        
        #get a dict[staxid][...] for each group
        if(len(E) > 0):
            self.E = get_top_bitscore_counts_per_taxid(E)
        if(len(K) > 0):
            self.K = get_top_bitscore_counts_per_taxid(K)
        if(len(R) > 0):
            self.R = get_top_bitscore_counts_per_taxid(R)
        if(len(U) > 0):
            self.U = get_top_bitscore_counts_per_taxid(U)
            
        #sort taxids on bitscore/count/numerical order
        self.txidsE = sorted(list(self.E.keys()), key=lambda txid: (self.E[txid]['bitscore'], self.E[txid]['count'], txid), reverse=True)
        self.txidsK = sorted(list(self.K.keys()), key=lambda txid: (self.K[txid]['bitscore'], self.K[txid]['count'], txid), reverse=True)
        self.txidsR = sorted(list(self.R.keys()), key=lambda txid: (self.R[txid]['bitscore'], self.R[txid]['count'], txid), reverse=True)
        self.txidsU = sorted(list(self.U.keys()), key=lambda txid: (self.U[txid]['bitscore'], self.U[txid]['count'], txid), reverse=True)
        #N doesn't need to be sorted because there will be only one or none
        self.txidsN = list(self.N.keys())
        
    #accessors
    def E(self):
        return(self.E)
    def K(self):
        return(self.K)
    def R(self):
        return(self.R)
    def U(self):
        return(self.U)
    def N(self):
        return(self.N)
        
    def txidsE(self):
        return(self.txidsE)
    def txidsK(self):
        return(self.txidsK)
    def txidsR(self):
        return(self.txidsR)
    def txidsU(self):
        return(self.txidsU)
    def txidsN(self):
        return(self.txidsN)
        
    def init_des(self):
        return(self.init_des)
    def set_final_des(self, final_des):
        self.final_des = final_des
    def final_des(self):
        return(self.final_des)
        
    def otu(self):
        return(self.otu)
    def size(self):
        return(self.size)
    def num_hits(self):
        return(self.num_hits)
        
    ### temp/debugging
    def set_choice(self, choice):
        self.choice = choice
    def choice(self):
        return(self.choice)
#

def get_top_bitscore_counts_per_taxid(lol):
        #   lol = [ [accnver, pident, evalue, bitscore, staxid, stitle, qcovs], [...] ]
        #returns: dict[staxid][bitscore|count|pident|stitle|accnver[]]
        
        #R[(141|100.000|3712|47) (141|100.000|3711|36) 
        #R(bitscore|pident|staxid|count)...
        
        d = {}
        staxids = []
        
        #for each staxid, get counts of each bitscore
        for hit in lol:
            #hit example:
            # ['LR031876.1', '97.368', '4.61e-27', '130', '3712', 'Brassica oleracea HDEM genome, scaffold: C7', 100]
            # [accnver,       pident,   evalue,  bitscore, staxid, stitle, qcovs]
            # [   0              1         2        3         4       5      6]
            #we don't need to keep 'evalue'
            staxid = hit[4]
            bitscore = float(hit[3])
            pident = float(hit[1])
            stitle = hit[5]
            accnver = hit[0]
            qcovs = hit[6]
            
            #initialize or update info for the current staxid
            if staxid not in d:
                d[staxid] = {}
                d[staxid]["bitscore"] = bitscore
                d[staxid]["count"] = 1
                d[staxid]["pident"] = pident #NOTE: the pident may not necessarily be the same for all identical staxid/bitscores pairs
                d[staxid]["stitle"] = stitle
                d[staxid]["accnver"] = [accnver]
                d[staxid]["qcovs"] = qcovs
            elif bitscore < d[staxid]["bitscore"]:
                #staxid is in d but bitscore is lower than the existing bitscore (because they have been sorted in descending order)
                pass
            else:
                #staxid is in d but bitscore is the same, so we count it, and append the current accnver
                d[staxid]["count"] += 1
                d[staxid]["accnver"] .append(accnver)
                
        return d
#

# class Error(Exception):
   # """Base class for other exceptions"""
   # pass
# class ExplicitException(Error):
   # """Raised when the input value is too small"""
   # pass

def parse_blast_hit(blasthit, indexOf, accnsWithDubiousTaxAssigns):
    '''Parse the list returned from the ParseBlastout iterator'''
    #TODO: replace multiple try/except blocks... pre-determine type elsewhere and send appropriate function with each call?
    #Query line (blasthit[0]) format examples:                                   #  otu ... size
    # Query: Otu1 1343563                                                        # Otu1 ... 1343563
    # Query: NB501124:253:HFJFVBGXB:1:11101:11597:1050 1:N:0:CACCGG;size=726951; # NB501124:253:HFJFVBGXB:1:11101:11597:1050 1:N:0:CACCGG; ... 726951
    # Query: TRINITY_DN393_c0_g1_i1 len=27077 path=[27055:0-27076]               # TRINITY_DN393_c0_g1_i1 ... 27077
    # Query: NODE_1_length_19883_cov_76.107931_g0_i0                             # NODE_1_length_19883_cov_76.107931_g0_i0 ... 19883
    failed = False
    try:
        #format 1 (usearch OTU + abundance)
        (otu, size) = blasthit[0][9:].split(' ')
    except:
        try:
            #format 2 (usearch-dereplicated fastq)
            (otu, size) = blasthit[0][9:].split('size=')
        except:
            try:
                #format 3 (trinity contig)
                bhit = re.sub(r' path=.+', '', blasthit[0][9:])
                (otu, size) = bhit.split(' len=')
            except:
                try:
                    #format 4 (Spades contig)
                    otu = blasthit[0][9:]
                    tmp = re.sub(r'_cov_.+', '', otu)
                    size = re.sub(r'NODE_\d+_length_', '', tmp)
                except:
                    failed = True
                    
    if failed:
        print("Error: Could not find an 'otu' and/or 'size' in this line:", file=sys.stderr)
        print(blasthit[0], file=sys.stderr)
        print("(From func: parse_blast_hit)", file=sys.stderr)
        sys.exit(0)
        
    if ';size' in size:
        size = re.sub(r';', '', size)
        size = size.split('=')[1]
    if ';' in size:
        size = re.sub(r';', '', size)
        
    maxbitscore = 0
    bitscores = []
    values = []
    lol = [] # lines of lines. temp var
    
    #Otu1_633460 [EU]->[K]
    #  E[(450|100.000|77133|144) (450|100.000|194843|1) (450|100.000|348578|1)
    
    #delete comment lines
    del blasthit[0:4]
    
    idx = 0
    for line in blasthit:
        hitvals = line.split('\t')
        # indexOf = {'stitle': 8, 'bitscore': 6, 'evalue': 5, 'mismatch': 4, 'pident': 2, 'length': 3, 'staxids': 7, 'sseqid': 1, 'qseqid': 0}
        #query id, subject id, % identity, alignment length, mismatches, evalue, bit score, subject tax ids, subject title
        # 0           1            2          3               4            5      6          7               8
        #['LR031899.1', '90.789', '76', '3', '1.30e-17', '99.0', '3712', 'Brassica oleracea HDEM' ]
        try:
            accnver = hitvals[indexOf['sseqid']].split('|')[-2]
        except:
            print("'accnver = hitvals[indexOf['sseqid']].split('|')[-2]'  failed on this line:")
            print(line)
            sys.exit(0)
        #skip this hit if it's accnver is on our list of dubiously classified sequences
        if accnver in accnsWithDubiousTaxAssigns:
            continue
            
        pident = hitvals[indexOf['pident']]
        evalue = hitvals[indexOf['evalue']]
        bitscore = hitvals[indexOf['bitscore']]
        staxids = hitvals[indexOf['staxids']]
        stitle = hitvals[indexOf['stitle']][:100]
        qcovs = hitvals[indexOf['qcovs']]
        
        #We need blast bit scores to be in descending order but they are frequently not (multi-part hits to the same subject),
        # so we create a list of lists here, for sorting later, so that downstream operations will not croak
        #append [bitscore, index] to bitscores
        bitscores.append( [bitscore, idx] )
        
        #get rid of OTU/seq name element
        hitvals = hitvals[1:]
        #store just what we want
        lol.append([accnver, pident, evalue, bitscore, staxids, stitle, qcovs])
        
        idx = idx + 1
        
    # reorder blasthit lines #
    #1) sort blasthit lines descending by bitscores
    bitscores.sort(key=lambda x: float(x[0]), reverse=True)
    #2) get the sorted indexes
    sortedIDXs = [bitscores[i][1] for i,j in enumerate(bitscores)]
    #3) resort lines
    lol = [lol[i] for i in sortedIDXs]
    
    if(0):
        if len(lol) > 1:
            print("OTU", otu)
            print("size", size)
            print('\n'.join(str(' '.join(str(x) for x in v)) for v in lol))
            sys.exit(0)
            
    return (otu, size, lol)
#


def get_Fields(opts):
    '''Reads up to the first "Fields" line, parses it, and returns a list'''
    # Fields: query id, subject id, % identity, alignment length, mismatches, evalue, bit score, subject seq
    with open(opts['-i']) as infile:
         #skip irrelevant lines
        for line in infile:
            if re.search("# Fields:", line) is None:
                continue
            else:
                break

        line = re.sub('# Fields: ', '', line.rstrip())
        line = re.sub('%', 'pct', line)
        line = re.sub(', ', '\t', line)
        line = re.sub(' ', '_', line)
        fields = line.split('\t')
        return fields
#


def parse_options(argv):
    #create an ArgumentParser object
    parser = argparse.ArgumentParser(
        description="Create lists of non-contaminating (target) and contaminating (non-target) OTUs by parsing local BLAST results of your sequences against the nt database.")
    
    #add arguments
    requiredArg = parser.add_argument_group('required arguments')
    requiredArg.add_argument("-i", "--input_fp", help="BLAST output file (using the format options given in #3 below)")
    requiredArg.add_argument("-k", "--keeper_txids", help="TaxonIDs of intended target organisms that should be KEPT.")
    requiredArg.add_argument("-e", "--environ_txids", help="TaxonIDs of Environmental samples of intended target organisms (provisionally KEPT. see below).")
    requiredArg.add_argument("-r", "--reject_txids", help="TaxonIDs of likely contaminant organisms that should be REJECTED.")
    requiredArg.add_argument("-m", "--merge_dmp", help="The merge.dmp file from NCBI to update taxIDs **Important if your local BLAST db is old but your taxIDs are new**")
    requiredArg.add_argument("-t", "--taxassndir", help="The directory that AccnsWithDubiousTaxAssigns.txt is or will be created in. Usually where USEARCH is located.")
    
    #parser.add_argument("-o", "--output_dir", help="The directory to save output files. Default: current directory")
    parser.add_argument("-p", "--prefix", help="Prefix to use for [otus2keep.txt] and [otus2filter.log]. Eg., 'CLas.CACCGG'")
    parser.add_argument("-n", "--num_hits", help="Number of hits to process (to terminate early)")
    parser.add_argument("-f", "--force_overwrite", help="Force overwrite of output files (including log).",
                        action='store_true')
    
    #parse the command line parameters
    opts = parser.parse_args()
    
    opts.dubiousIDs = os.path.join(opts.taxassndir, "AccnsWithDubiousTaxAssigns.txt")
    
    return opts
#


def all_input_files_exist(opts):
    '''Checks whether all input files exist. Returns: bool'''
    is_error = False
    not_exist = []

    #blastout file
    if not os.path.isfile(opts.input_fp):
        print ("** Error ** Input file ["+opts.input_fp+"] does not exist!")
        is_error = True
    #taxID files
    for f in opts.keeper_txids.split(","):
        if not os.path.isfile(f):
            print ("** Error ** Input file ["+f+"] does not exist!")
            is_error = True
    for f in opts.environ_txids.split(","):
        if not os.path.isfile(f):
            print ("** Error ** Input file ["+f+"] does not exist!")
            is_error = True
    for f in opts.reject_txids.split(","):
        if not os.path.isfile(f):
            print ("** Error ** Input file ["+f+"] does not exist!")
            is_error = True
    #merge.dmp file
    if not os.path.isfile(opts.merge_dmp):
        print ("** Error ** Input file ["+opts.merge_dmp+"] does not exist!")
        is_error = True
    #taxassndir directory
    if not os.path.isdir(opts.taxassndir):
        print ("** Error ** Input directory ["+opts.taxassndir+"] does not exist!")
        is_error = True

    if is_error:
        return False
    else:
        return True
#


def load_IDs(opts):
    '''Loads taxonIDs of K, E and R groupings. Returns: dict'''
    ids_d = {}
    for f in opts.keeper_txids.split(","):
        print("Loading K IDs in ["+f+"]")
        with open(f, 'r') as idfile:
            for line in idfile:
                taxid = int(line.rstrip())
                ids_d[taxid] = "K"

    for f in opts.environ_txids.split(","):
        print("Loading E IDs in ["+f+"]")
        with open(f, 'r') as idfile:
            for line in idfile:
                taxid = int(line.rstrip())
                ids_d[taxid] = "E"

    for f in opts.reject_txids.split(","):
        print("Loading R IDs in ["+f+"]")
        with open(f, 'r') as idfile:
            for line in idfile:
                taxid = int(line.rstrip())
                ids_d[taxid] = "R"

    #phiX
    ids_d[10847] = "R"

    return ids_d
#


def load_accnsWithDubiousTaxAssigns(opts):
    print("Loading [AccnsWithDubiousTaxAssigns.txt]")
    d = {}
    with open(opts.dubiousIDs, 'r') as file:
        for line in file:
            if '\t' in line:
                n = line.rstrip().split('\t')[0]
                d[n] = 1
    return(d)
#


def add_old_txids_to_EKRs(usersTaxIDs, mergedDmpIDs):
    '''Adds "old" taxonomic IDs to the (newer) userTaxIDs to ensure any old blastout IDs are classified into E/K/Rs correctly'''
    #mergedDmpIDs[newID] = oldID
    #usersTaxIDs[newID] = E|K|R
    
    #usersTaxIDs ARE (should be!) NEWER THAN THE BLAST DB
    #SO WE NEED TO KNOW WHAT TO DO with THE OLD BLAST taxIDs
    #We need to UPDATE usersTaxIDs TO CONTAIN OLD IDs
    num_oldIDs_added = 0
    for newID in usersTaxIDs.keys():
        if newID in mergedDmpIDs:
            #THEN THERE IS AN oldID ALSO
            oldID = mergedDmpIDs[newID]
            #so we add the OLD ID and set it to the SAME EKR GROUPING AS THE NEW ID
            usersTaxIDs[oldID] = usersTaxIDs[newID] # E|K|R
            num_oldIDs_added += 1
            
    print("[" + str(num_oldIDs_added) + "] newIDs had an oldID in merge.dmp and were added to the E/K/R usersTaxID lists", file=sys.stderr)
    
    return usersTaxIDs
#


def load_mergedDmpIDs(opts):
    print("Loading [merge.dmp]")
    mergedDmpIDs = {}
    with open(opts.merge_dmp, 'r') as file:
        for line in file:
            splitline = line.split('\t') #281694	|	369572	|
            (oldID, newID) = [ splitline[i] for i in (0,2) ]
            mergedDmpIDs[newID] = oldID
            
    return mergedDmpIDs
#


def blastHeaders():
    #@headrIDX{@blastheaders} = @blastterms;
    blastheaders = ["query id","subject id","% identity","alignment length","mismatches","evalue","bit score","subject tax ids","subject title","% query coverage per subject"]
    blastterms = ["qseqid","sseqid","pident","length","mismatch","evalue","bitscore","staxids","stitle","qcovs"]
    headrIDX = {}
    for i,orighdr in enumerate(blastheaders):
        headrIDX[orighdr] = blastterms[i]
    
    return headrIDX
#


def getIndexOf(opts, headrIDX):
    '''Get indexes of required headers in the the current blastout file. Returns dict'''
    indexOf = {}
    with open(opts.input_fp, 'r') as file:
        for line in file:
            if '# Fields: ' in line:
                line = line[10:]
                blastheaders = line.rstrip().split(', ')
                #convert blast's 'Fields:' output line into blastn's format versions for the command line (eg. '% identity' -> 'pident')
                for i,orighdr in enumerate(blastheaders):
                    if orighdr in headrIDX:
                        indexOf[headrIDX[orighdr]] = i
                break
                
        return indexOf
#


def keep_or_reject_judgement1(bh, scrdes):
    '''Use info from the blast hits to decide whether an OTU/sequence should be [K]ept or [R]ejected'''
    
    #concatenate whatever E/K/R/U/N groups that were found in this blast hit - the "initial designation"
    init_des = ''.join(bh.init_des)
    
    #TODO: just use bh.txidsE directly?
    txidsE = bh.txidsE
    txidsK = bh.txidsK
    txidsR = bh.txidsR
    txidsU = bh.txidsU
    
    #N doesn't need to be sorted because there will be only one or none
    txidsN = list(bh.N.keys())
    
    #find the best bitscore and its number of occurences within in each group
    #TODO: calculate/use total unique accns at the highest bitscore? This would be > len(txidsE)
    if len(txidsE) > 0:
        highestbitE = bh.E[txidsE[0]]['bitscore']
        totalhighestbitcountE = 0
        numberOfHighestBitScoresE = 0
        for txid in txidsE:
            if bh.E[txid]['bitscore'] == highestbitE:
                totalhighestbitcountE += bh.E[txid]['count']
                numberOfHighestBitScoresE += 1
                
    if len(txidsK) > 0:
        highestbitK = bh.K[txidsK[0]]['bitscore']
        totalhighestbitcountK = 0
        numberOfHighestBitScoresK = 0
        for txid in txidsK:
            if bh.K[txid]['bitscore'] == highestbitK:
                totalhighestbitcountK += bh.K[txid]['count']
                numberOfHighestBitScoresK += 1
                
    if len(txidsR) > 0:
        highestbitR = bh.R[txidsR[0]]['bitscore']
        totalhighestbitcountR = 0
        numberOfHighestBitScoresR = 0
        for txid in txidsR:
            if bh.R[txid]['bitscore'] == highestbitR:
                totalhighestbitcountR += bh.R[txid]['count']
                numberOfHighestBitScoresR += 1
                
    if len(txidsU) > 0:
        highestbitU = bh.U[txidsU[0]]['bitscore']
        totalhighestbitcountU = 0
        numberOfHighestBitScoresU = 0
        for txid in txidsU:
            if bh.U[txid]['bitscore'] == highestbitU:
                totalhighestbitcountU += bh.U[txid]['count']
                numberOfHighestBitScoresU += 1
                
                
                
    #look for misclassified accession numbers in K or R when highest bitscores are equal
    if ("K" in init_des and "R" in init_des) and (highestbitK == highestbitR):
        if totalhighestbitcountK == 1 and totalhighestbitcountR > 20:
            accnK = bh.K[txidsK[0]]['accnver'][0]
            if accnK not in singletKaccns:
                singletKaccns[accnK] = {}
                singletKaccns[accnK]['count'] = 1
                singletKaccns[accnK]['taxid'] = txidsK[0]
            else:
                singletKaccns[accnK]['count'] += 1
        elif totalhighestbitcountR == 1 and totalhighestbitcountK > 20:
            accnR = bh.R[txidsR[0]]['accnver'][0]
            if accnR not in singletRaccns:
                singletRaccns[accnR] = {}
                singletRaccns[accnR]['count'] = 1
                singletRaccns[accnR]['taxid'] = txidsR[0]
            else:
                singletRaccns[accnR]['count'] += 1
                
    ###
    # Equivalent meaning/values between Perl and Python code
    # Perl                           Python                     # Comment...
    # $des{E}{$taxidsE[0]}{C}        bh.E[txidsE[0]]['count']   # number of hits for the top taxid (but we're using totalhighestbitcountE, which seems more logical)
    # $des_nHits{E}                  len(txidsE)                # number of unique taxid hits
    # $txID_bestBit{$taxidsE[0]}{B}  highestbitE                # the highest bitscore in E
    # NA                             totalhighestbitcountE      # number of hits having the highest bitscore (across all taxids with best bitscore)
    ###
    #make final Keep or Reject decision
    choice = []
    if init_des == 'EKRU': #Desc1
        final_des = "K" #default
        choice = [0]
        if highestbitK < highestbitR:
            final_des = "R"
            choice = [1, len(txidsK), len(txidsR), totalhighestbitcountK, totalhighestbitcountR, "|", numberOfHighestBitScoresK, numberOfHighestBitScoresR, singletRaccns]
        #elif highestbitK == highestbitR and bh.K[txidsK[0]]['count'] <= bh.R[txidsR[0]]['count']: ##<- this is the same as the Perl code but we've decided not to use it
        elif highestbitK == highestbitR and totalhighestbitcountK <= totalhighestbitcountR and len(txidsK) <= len(txidsR): ## <- experimental... but we think better
            final_des = "R"
            choice = [2, len(txidsK), len(txidsR), totalhighestbitcountK, totalhighestbitcountR, "|", numberOfHighestBitScoresK, numberOfHighestBitScoresR, singletRaccns]
            #Otu32751_26 [EKRU]->[R] #<- very interesting... check w/james
            #Otu54198_14 [EKRU]->[R] #<- also interesting... check that the highest pidents are being kept/not overwritten with lower ones
            
        elif highestbitK == highestbitU and totalhighestbitcountK <= totalhighestbitcountU: ##
            final_des = "R"
            choice = [3, len(txidsK), len(txidsR), totalhighestbitcountK, totalhighestbitcountR, "|", numberOfHighestBitScoresK, numberOfHighestBitScoresR, singletRaccns]
            
    elif init_des == 'EKR':  #Dec2
        final_des = 'K'
        if highestbitK < highestbitR:
            final_des = "R"
            choice = [1, highestbitK, highestbitR, len(txidsK), len(txidsR), singletRaccns]
        elif highestbitK == highestbitR and totalhighestbitcountK <= totalhighestbitcountR: ##
            final_des = "R"
            choice = [2, len(txidsK), len(txidsR), totalhighestbitcountK, totalhighestbitcountR, singletRaccns]
            
    elif init_des == 'EU':   #Dec3
        final_des = 'K'
        # elsif ($txID_bestBit{$taxidsE[0]}{B} < $txID_bestBit{$taxidsU[0]}{B}  && $des{E}{$taxidsE[0]}{C} < $des{U}{$taxidsU[0]}{C})
        if highestbitE < highestbitU and totalhighestbitcountE <= totalhighestbitcountU:
            final_des = "R"
            choice = [1, len(txidsE), len(txidsU), totalhighestbitcountE, totalhighestbitcountU]
        # elsif ($txID_bestBit{$taxidsE[0]}{B} < $txID_bestBit{$taxidsU[0]}{B}  #if(best_bitscore for E < best_bitscore U
            # && $des_nHits{E} < $des_nHits{U})                                 # AND if num_hits_for E < num_hits_for U)
        elif highestbitE < highestbitU and len(txidsE) < len(txidsU):
            final_des = "R"
            choice = [2, len(txidsE), len(txidsU), totalhighestbitcountE, totalhighestbitcountU]
        # elsif ($txID_bestBit{$taxidsE[0]}{B} == $txID_bestBit{$taxidsU[0]}{B} && $stitle[3] =~ /[Cc]hloroplast/)   #OR if(stitle[3] =~ /[Cc]hloroplast/)
        elif highestbitE == highestbitU and bool(re.search('[Cc]hloroplast', bh.U[txidsU[0]]['stitle'])):
            final_des = "R"
            choice = [3, len(txidsE), len(txidsU), totalhighestbitcountE, totalhighestbitcountU]
        
    elif init_des == 'KU':   #Dec4
        final_des = 'K'
        # elsif ($txID_bestBit{$taxidsK[0]}{B} <  $txID_bestBit{$taxidsU[0]}{B}  #if(best_bitscore for K < best_bitscore U
            # && $des_nHits{K} < $des_nHits{U})                                 # AND if num_hits_for K < num_hits_for U)
        if highestbitK < highestbitU and len(txidsK) < len(txidsU):
            final_des = "R"
            choice = [1, len(txidsK), len(txidsU)]
        # elsif ($txID_bestBit{$taxidsK[0]}{B} == $txID_bestBit{$taxidsU[0]}{B} && $des_nHits{K} <= $des_nHits{U})
        elif highestbitK == highestbitU and len(txidsK) <= len(txidsU):
            final_des = "R"
            choice = [2, len(txidsK), len(txidsU)]
        
    elif init_des == "KRU":  #Dec5
        final_des = 'K'
        # choice = [0, phiX_txid, txidsU]
        if highestbitK < highestbitR:
            final_des = "R"
            choice = [1, len(txidsK), len(txidsR), len(txidsU), singletRaccns]
        elif highestbitK <= highestbitU and phiX_txid in txidsR:
            final_des = "R"
            choice = [2, len(txidsK), len(txidsR), len(txidsU), totalhighestbitcountK, totalhighestbitcountR, totalhighestbitcountU, singletRaccns]
        elif highestbitK == highestbitR and totalhighestbitcountK <= 2 and totalhighestbitcountK < totalhighestbitcountR:
            final_des = "R"
            choice = [3, len(txidsK), len(txidsR), len(txidsU), totalhighestbitcountK, totalhighestbitcountR, totalhighestbitcountU, singletRaccns]
        elif highestbitK == highestbitR or (phiX_txid in txidsR and totalhighestbitcountR > 2):
            final_des = "R"
            choice = [4, len(txidsK), len(txidsR), len(txidsU), totalhighestbitcountK, totalhighestbitcountR, totalhighestbitcountU, singletRaccns]
        
    elif init_des == "KR":  #Dec6
        final_des = 'K'
        if highestbitK < highestbitR:
            final_des = "R"
            choice = [1, highestbitK, highestbitR, len(txidsK), len(txidsR), singletRaccns]
        elif highestbitK == highestbitR and totalhighestbitcountK <= totalhighestbitcountR:
            final_des = "R"
            choice = [2, len(txidsK), len(txidsR), totalhighestbitcountK, totalhighestbitcountR, singletRaccns]
            
    else:
        if init_des in scrdes:
            final_des = scrdes[init_des]
        else:
            print("*** Error *** Designation", init_des, "not found in designator dictionary 'scrdes'")
            sys.exit(0)
            
    bh.set_final_des(final_des)
    bh.set_choice(choice)
#


def print_to_log(bh, logfh):
    #print OTU name, size and initial designation
    logfh.write('_'.join([bh.otu, bh.size]) + ' ['+''.join(bh.init_des)+']->' + '['+bh.final_des+']\n')
    
    #print results
    logfh.write('  E[')
    if len(bh.txidsE) > 0:
        logfh.write(formated_hits(bh.E, bh.txidsE))
    else:
        logfh.write(']\n')
        
    logfh.write('  K[')
    if len(bh.txidsK) > 0:
        logfh.write(formated_hits(bh.K, bh.txidsK))
    else:
        logfh.write(']\n')
        
    logfh.write('  R[')
    if len(bh.txidsR) > 0:
        logfh.write(formated_hits(bh.R, bh.txidsR))
    else:
        logfh.write(']\n')
        
    logfh.write('  U[')
    if len(bh.txidsU) > 0:
        logfh.write(formated_hits(bh.U, bh.txidsU))
    else:
        logfh.write(']\n')
        
    if len(bh.txidsN) > 0:
        logfh.write('  N[(')
        logfh.write('|'.join(['', '', bh.txidsN[0],'']) + ") ]\n")
        
    # #print choice
    # if len(bh.choice) > 0:
        # logfh.write('  '+''.join(bh.init_des)+'choice: '+ str(bh.choice) + '\n')
        
    logfh.write('\n')
#


def tabulate_results_summary(bh, des_summary):
    #tabulate results summary
    init_des = ''.join(bh.init_des)
    final_des = ''.join(bh.final_des)
    des_summary[init_des][final_des]['OTUs'] += 1
    des_summary[init_des][final_des]['Reads'] += int(bh.size)
#


def save_designation_summary(outf4, des_summary):
    with open(outf4, "w") as sumfh:
        d = json.loads(json.dumps(des_summary))
        initial_designators = sorted(d.keys())
        sumfh.write(('{:^7}|{:^8}|{:^12}|{:^10}|{:^10}'.format("Dsgntr","OTUs","Reads","Rd/OTU","newDsgntr")))
        sumfh.write("\n")
        tot_num_otusK = 0
        tot_num_otusR = 0
        tot_num_otus = 0
        tot_num_readsK = 0
        tot_num_readsR = 0
        tot_num_reads = 0
        for init_des in initial_designators:
            finaldes = sorted(d[init_des].keys())
            for final_des in finaldes:
                num_otus  = d[init_des][final_des]['OTUs']
                num_reads = d[init_des][final_des]['Reads']
                reads_per_otu = int(round(num_reads / num_otus))
                sumfh.write(('{:<7}|{:>7} |{:>11} |{: >9} |{: ^10}'.format(init_des, num_otus, num_reads, reads_per_otu, final_des)))
                sumfh.write("\n")
                #calc totals
                tot_num_otus  += num_otus
                tot_num_reads += num_reads
                #calc totals by final designation
                if final_des == 'K':
                    tot_num_otusK  += num_otus
                    tot_num_readsK += num_reads
                else:
                    tot_num_otusR  += num_otus
                    tot_num_readsR += num_reads
        #
        sumfh.write(('{:<7}|{:>7} |{:>11} |{: >9} |{: ^10}'.format("TotalK:", tot_num_otusK, tot_num_readsK, "", "K")))
        sumfh.write("\n")
        sumfh.write(('{:<7}|{:>7} |{:>11} |{: >9} |{: ^10}'.format("TotalR:", tot_num_otusR, tot_num_readsR, "", "R")))
        sumfh.write("\n")
        sumfh.write(('{:<7}|{:>7} |{:>11} |{: >9} |{: ^10}'.format("Total:", tot_num_otus, tot_num_reads, "", "K+R")))
        sumfh.write("\n")
        
#


def save_dubiously_assigned_accessions(outf5, singletRaccns, singletKaccns):
    #sort accns decending on abundance
    accnRs = sorted(list(singletRaccns.items()), key=lambda x: x[1]['count'], reverse=True)
    accnKs = sorted(list(singletKaccns.items()), key=lambda x: x[1]['count'], reverse=True)
    
    #save to file
    with open(outf5, "w") as badfh:
        badfh.write("#Bad[R]accns\tCount\tTaxID\n")
        for txid, d in accnRs:
            badfh.write(txid + "\t" + str(d['count']) + "\t" + str(d['taxid']) + "\n")
        badfh.write("\n#Bad[K]accns\tCount\tTaxID\n")
        for txid, d in accnKs:
            badfh.write(txid + "\t" + str(d['count']) + "\t" + str(d['taxid']) + "\n")
#


def formated_hits(bhX, txidsX):
    #return a blank if there's nothing to print
    if len(txidsX) == 0:
        return ''
        
    #concatenate hit info for all hits into a single line
    hit = ''
    for txid in txidsX:
        bitscore = '{0:g}'.format(bhX[txid]['bitscore'])
        hit = hit + '(' + '|'.join([bitscore, str(bhX[txid]['pident']), str(txid), str(bhX[txid]['count']), str(bhX[txid]['qcovs']), ]) + ') '
        
    #concatenate just the subject title on the next line
    hit = hit + '\n   (' + bhX[txidsX[0]]['stitle'] + ') ]\n'
    
    return hit
#


def main(argv):
    '''main function'''
    opts = parse_options(argv)
    
    if not all_input_files_exist(opts):
        sys.exit(0)
        
    des_summary = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    accnsWithDubiousTaxAssigns = load_accnsWithDubiousTaxAssigns(opts)
    mergedDmpIDs = load_mergedDmpIDs(opts)
    usersTaxIDs_d = load_IDs(opts)
    usersTaxIDs_d = add_old_txids_to_EKRs(usersTaxIDs_d, mergedDmpIDs)
    
    #create a dictionary of classification decisions. values with "Dec1", etc. require more complex decision-making code
    scrdes = dict(list(zip("E EK EKU K KU  ER ERU KR KRU R RU U  EKRU EKR EU".split(), "K K K K Dec4  R R Dec6 Dec5 R R R  Dec1 Dec2 Dec3".split())))
    
    if(1): #TODO: make this a default option in opts
        scrdes["N"] = "K"
    else:
        scrdes["N"] = "R"
        
    headrIDX = blastHeaders()
    indexOf = getIndexOf(opts, headrIDX)
    
    #create a parser object
    parser = ParseBlastout(opts.input_fp)
    
    #prep output filenames
    outf1 = "otus2filter.log"
    outf2 = "otus2keep.txt"
    outf3 = "otus2reject.txt"
    outf4 = "otus2summary.txt"
    outf5 = "bad_accns.txt"
    if opts.prefix:
        if opts.prefix[-1] != ".":
            opts.prefix = opts.prefix + "."
        outf1 = opts.prefix + outf1
        outf2 = opts.prefix + outf2
        outf3 = opts.prefix + outf3
        outf4 = opts.prefix + outf4
        outf5 = opts.prefix + outf5
        
    #open output files
    with open(outf1, "w") as logfh, open(outf2, "w") as kfh, open(outf3, "w") as rfh:
        
        #initialize a dict to check whether we encounter the same OTU more than once
        otus = {}
        #parse the blastout file
        num_hits = 0
        for blasthit in parser:
            #parse the blasthit object
            (otu, size, hitlines) = parse_blast_hit(blasthit, indexOf, accnsWithDubiousTaxAssigns)
            
            if otu in otus:
                #we've got a concatenated blastout file with re-blasted OTUs
                #TODO add/use this info instead of skipping?
                # unecessary IF we put the deeper/better blasts at the top. if not, we're tossing good info
                # but for now, we'll skip
                continue
            otus[otu] = 1
            
            #create a blastHit object
            bh = blastHit(otu, size, hitlines, usersTaxIDs_d)
            
            #make classification decision
            keep_or_reject_judgement1(bh, scrdes)
            
            #add results to summary
            tabulate_results_summary(bh, des_summary)
            
            #save keeper/reject IDs to file
            if bh.final_des == "K":
                kfh.write(otu + "\n")
            elif bh.final_des == "R":
                rfh.write(otu + "\n")
                
            #print results to log as we go
            print_to_log(bh, logfh)
            
            num_hits += 1
            if opts.num_hits and num_hits == int(opts.num_hits):
                break
                
    #save other tabulated info
    save_designation_summary(outf4, des_summary)
    save_dubiously_assigned_accessions(outf5, singletRaccns, singletKaccns)
#


if __name__ == "__main__":
    main(sys.argv)