#!/usr/bin/env python

"""
compare_categories2matrix.py

Parses pairwise adonis/anosim output files into matrices for easier viewing

Usage:
    compare_categories2matrix.py -m <metric> [-t <type>] [-d <dir>] [-o <filepath>] [-r <ret_type>]
    compare_categories2matrix.py --version
    compare_categories2matrix.py -h

Arguments:
    -m <metric>         Metric used (eg. 'hellinger'). See note below

Options:
    -d <dir>            Directory where the anosim/adonis files are located [default: ./]
    -t <type>           Stat type ('anosim' or 'adonis'). Auto-determined if not defined
    -o <filepath>       Output filepath. Prints to stdout if not defined
    -r <ret_type>       Return value type ('pvalue' or 'R2'). [default: pvalue]
    -h --help           This help
    --version           Version

Note: input adonis/anosim filenames MUST be in this form:
<metric>_<category>__<pair1>_v_<pair2>_results.txt
Eg.:
hellinger_TissueType__Budwood_v_Leaf_results.txt

Related scripts (customized versions):
1) compare_categories.py (used to create pairwise anosim results files)
2) adonis.r              (used to create pairwise adonis results files)
"""


from __future__ import print_function
from os.path import isfile
from docopt import docopt
import sys
import os
import re
import pandas as pd
from StringIO import StringIO



def get_and_check_opts(args):
    #docopt(doc, argv=None, help=True, version=None, options_first=False)
    opts = docopt(__doc__, version="compare_categories2matrix.py, version 0.3")
    
    #check opts
    message = ''
    optionsFail = False
    if not opts['-m']:
        optionsFail = True
        message = "*** Error *** -m <metric> not defined\n"
    
    #ensure dir ends with a '/' char
    if opts['-d'][-1] != '/':
        opts['-d'] = opts['-d'] + '/'
    
    #if -t is given, make sure its on the list
    if opts['-t'] in ['anosim', 'adonis']:
        pass
    elif opts['-t'] is not None:
        optionsFail = True
        message = (message+"*** Unknown stat type ["+opts['-t']+"] ***\n"
                    "Known types are 'anosim' and 'adonis'\n")
    
    #if -r is given, make sure its on the list
    if opts['-r'] in ['pvalue', 'R2']:
        pass
    else:
        optionsFail = True
        message = (message+"*** Unknown return type requested ["+opts['-r']+"] ***\n"
                    "Known types are 'pvalue' and 'R2'\n")
    
    if optionsFail:
        print(message)
        sys.exit(0)
        
    #add an underscore to metric if not provided (for regex)
    if opts["-m"][-1] == "_":
        opts["-m2"] = opts["-m"][:-1]
    else:
        opts["-m2"] = opts["-m"]
        opts["-m"] = opts["-m"]+"_"
        
    return opts


def get_filenames(opts):
    #get files starting with -m's value
    files = []
    for file in os.listdir(opts["-d"]):
        if re.search(r"^"+opts["-m"], file) is not None:
            files.append(opts["-d"]+file)
            
    if len(files) == 0:
        print("*** Error *** No files staring with "+opts["-m"]+" could be found in ["+opts["-d"]+"]")
        sys.exit(0)
        
    return files


def parse_filename(opts, filename):
    #remove directory, if any
    filename = filename.split(os.sep)[-1]
    #remove metric from filename
    fname = re.sub(opts["-m"], '', filename)
    #remove end of name
    fname = re.sub(opts['--tail'], '', fname)
    #split what's left
    categ,pair = re.split("__", fname)
    #get pairs
    pair1,pair2 = re.split("_v_", pair)
    return [categ,pair1,pair2]


def get_initialized_dict_of_categories_from_filenames(opts, filenames):
    #eg: hellinger_Location__Indiantown_v_Weirsdale_results.txt
    categs_dict = {}
    for filename in filenames:
        categ,pair1,pair2 = parse_filename(opts, filename)
        
        #store in a dict#
        #initialize categ key (level 1)
        if categ not in categs_dict:
            categs_dict[categ] = {}
        #initialize pair1 and pair2 keys (level 2)
        if pair1 not in categs_dict[categ]:
            categs_dict[categ][pair1] = {}
        if pair2 not in categs_dict[categ]:
            categs_dict[categ][pair2] = {}
        if pair1 not in categs_dict[categ][pair2]:
            #initialize pair1 (level 3)
            categs_dict[categ][pair2][pair1] = 1
        if pair2 not in categs_dict[categ][pair1]:
            #initialize pair2 (level 3)
            categs_dict[categ][pair1][pair2] = 1
        #init selfs
        categs_dict[categ][pair1][pair1] = 0
        categs_dict[categ][pair2][pair2] = 0
            
    return categs_dict


def update_dict_from_anosim_result_files(opts, filenames, categs_dict):
    #Example input file:#
    # method name	ANOSIM
    # test statistic name	R
    # sample size	60
    # number of groups	2
    # test statistic	0.3754329501915708
    # p-value	0.001
    # number of permutations	999
    
    for filename in filenames:
        categ,pair1,pair2 = parse_filename(opts, filename)
        
        with open(filename, 'r') as infile:
            for line in infile:
                if re.search("^p-value", line) is not None:
                    x,pvalue = line.rstrip().split("\t")
                    break
                    #continue
                    
        #add pvalue to dict in both possible orders
        categs_dict[categ][pair1][pair2] = pvalue
        categs_dict[categ][pair2][pair1] = pvalue
        
    return categs_dict


def update_dict_from_adonis_result_files(opts, filenames, categs_dict):
    #Example input file:#
    # Call:
    # adonis(formula = as.dist(qiime.data$distmat) ~ qiime.data$m[[opts$category]],      permutations = opts$num_permutations) 

    # Permutation: free
    # Number of permutations: 999

    # Terms added sequentially (first to last)

    #                               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    # qiime.data$m[[opts$category]]  1    0.1002 0.100200  1.6464 0.01652  0.201
    # Residuals                     98    5.9641 0.060859         0.98348       
    # Total                         99    6.0643                  1.00000       
    for filename in filenames:
        categ,pair1,pair2 = parse_filename(opts, filename)
        
        with open(filename, 'r') as infile:
            for line in infile:
                if re.search("^qiime", line) is not None:
                    vals = re.split('\s+', line.rstrip())
                    if len(vals) >= 7:
                        pvalue = vals[6]
                        r2value = vals[5]
                        break
                    else:
                        #the qiime..[[opts$category]] line was split in two so we look for the next one
                        ### 200203. i can't find/create an adonis output file where the 'qiime.data$m...' line is continued
                        ###  on a second line, so i don't know how to handle 'R2' output. this code detects case and stops.
                        if opts['-r'] in 'R2':
                            print("** Sorry. I'm not sure I can find the 'R2' value in this multi-qiime.data-line adonis file **")
                            print("Please examine it's contents and modify [update_dict_from_adonis_result_files] as needed")
                            print("See file ["+filename+"]")
                            sys.exit(0)
                        ###
                        for line in infile:
                            if re.search("^qiime", line) is not None:
                                vals = re.split('\s+', line.rstrip())
                                if len(vals) >= 2:
                                    pvalue = vals[-1]
                                    #NOTE: some 'R2' code probably goes here IF the R2 value is found in this line
                                    break
                                else:
                                    print(line)
                                    print("** Oops ** function [update_dict_from_adonis_result_files] couldn't handle the above line")
                                    print("See file ["+filename+"]")
                                    sys.exit(0)
        #simply overwrite pvalue with r2value if user has requested R^2 values
        if opts['-r'] in 'R2':
            pvalue = r2value
        #add pvalue to dict in both possible orders
        categs_dict[categ][pair1][pair2] = pvalue
        categs_dict[categ][pair2][pair1] = pvalue
        
    return categs_dict


def make_matrix(opts, categs_dict, categ):
    #create a matrix of values from a given category
    names = sorted(categs_dict[categ].keys())
    #create an empty dataframe with named rows and columns
    df = pd.DataFrame(columns=names, index=names)
    for name in names:
        df.loc[name] = categs_dict[categ][name]
        
    return(df)


def print_to_stdout(opts, categs_dict):
    #print to stdout by way of StringIO to keep tabs in output
    orig_stdout = sys.stdout
    stringIOout = StringIO()
    sys.stdout = stringIOout
    for categ in sorted(categs_dict.keys()):
        df = make_matrix(opts, categs_dict, categ)
        print(categ)
        df.to_csv(stringIOout, sep="\t")
        print('')
        
    sys.stdout = orig_stdout
    print(opts['-t'] + " (" + opts['-r'] + ")")
    print(opts['-m2'])
    print(stringIOout.getvalue(), end="")


def print_to_file(opts, categs_dict):
    print("print_to_file not yet implemented")
    sys.exit(0)
    #TODO: make the copy-pasted code below work
    # with open(filename, 'w') as outfile:
        # for categ in sorted(categs_dict.keys()):
            # df = make_matrix(opts, categs_dict, categ)
            # outfile.write(categ)
            # df.to_csv(stringIOout, sep="\t")
            # print('')
            
        # print(opts['-t'])
        # print(opts['-m2'])
        # print(stringIOout.getvalue(), end="")
        # outfile.write("#OTUID\ttaxonomy\tconf\n")
        # outfile.write("\n".join(output))
    # print("File saved as ["+opts.output_fp+"]")


def determine_stat_type(opts, filename):
    stat_type = None
    with open(filename, 'r') as infile:
        for line in infile:
            if re.search("ANOSIM", line) is not None:
                stat_type = 'anosim'
                break
            elif re.search("^adonis\(", line) is not None:
                stat_type = 'adonis'
                break
                
    if stat_type is None:
        print("*** Error *** Unrecognized format: ["+filename+"]")
        sys.exit(0)
        
    return stat_type


def main(args):
    opts = get_and_check_opts(args)
    opts['--tail'] = '_results.txt'
    
    filenames = get_filenames(opts)
    categs_dict = get_initialized_dict_of_categories_from_filenames(opts, filenames)
    
    if opts['-t'] is None:
        opts['-t'] = determine_stat_type(opts, filenames[0])
    
    if opts['-t'] == 'anosim':
        categs_dict = update_dict_from_anosim_result_files(opts, filenames, categs_dict)
    elif opts['-t'] == 'adonis':
        categs_dict = update_dict_from_adonis_result_files(opts, filenames, categs_dict)
        
    if opts['-o'] is None:
        print_to_stdout(opts, categs_dict)
    else:
        print_to_file(opts, categs_dict)


if __name__ == "__main__":
    main(__doc__)