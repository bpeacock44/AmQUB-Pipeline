#!/usr/bin/env python

"""
This script takes a biom-formatted ASV table and multiplies or divides all values by 
a single user-defined value, or adjusts each ASV separately by importing a file
containing the amount of adjustment desired for each ASV.

Note: Addition and subtraction are not currently supported
"""

from __future__ import print_function, division
import sys, os
from os.path import isfile
import numpy as np
from biom import Table, load_table
from qiime.util import make_option, write_biom_table, parse_command_line_parameters


SCRIPT_INFO = {}
SCRIPT_INFO['brief_description'] = 'Perform multiplication, division or normalization on an ASV table in biom format'
SCRIPT_INFO['script_description'] = """biom_table_math_ops.py

This script can perform three types of adjustments on the values of a biom-formatted ASV table:
1) Multiply (-m) or divide (-d) all values of all ASVs by a SINGLE, user-supplied value.
2) Divide all values of each ASV by a LIST of user-supplied values: one for each ASV.
3) Normalize the table to unity (columns/samples will sum to 1).

Notes:
- Addition and subtraction are not supported.
- For adjustments of type 1), if an output file (-o) is not supplied, the script will only 
  DISPLAY the 'MAXIMUM POSSIBLE SUBSAMPLING DEPTH' (MPSD) of the input biom file, then it exits 
  without saving an output file. This allows one to quickly determine the MPSD for a biom file 
  without creating a new file (though the MPSD is also displayed after creating a new biom file).
- To perform adjustments of type 2), an additional input file is required that contains a list of
  ASV names and their corresponding adjustment values. This input file can be created with the
  [rrna_copy_adjustments.py] script.
- Note that the MPSD is typically less than the minimum sum of all sample sums *if the table has 
  fractional values*. It is calculated by using the same math functions employed in certain QIIME 
  scripts that cast values to integers rather than leaving them as floating-point numbers. Knowing
  the MPSD is useful when, for instance, performing subsampling in QIIME on a table containing 
  fractional values. Using a value <= to the MPSD will ensure that no samples are discarded. See 
  the '--depth' option in QIIME's [multiple_rarefactions_even_depth.py] script for more info. It
  is also worth noting that before performing multiple rarefactions, or a single subsampling, one 
  should multiply the ASV table by a large constant (eg., 100 or 1000) so that when QIIME truncates
  fractions it will only be a small portion of each value. Then, after rarefaction/subsampling, 
  undo the multiplication by dividing the resulting ASV table by the amount used for multiplication.

Related scripts:
  convert_RDP_taxa_to_gg_format.py
  rrna_copy_adjustments.py
  biom_table_math_ops.py"""


SCRIPT_INFO['required_options'] = [
    make_option('-i', '--input_fp', type="existing_filepath",
                help='The input biom ASV table filepath. Used only with --multiply_by option '),
]
SCRIPT_INFO['optional_options'] = [
    make_option('-o', '--output_fp', type="new_filepath", help='Path to store result file'),
    make_option('-m', '--multiply_by', type='float',
                help='Multiply ASV table values by this amount. Note: intended ' +
                'for use on a single ASV table, prior to rarefactions, to overcome rounding ' +
                'errors introduced by [multiple_rarefactions_even_depth.py] if operating on ' +
                'an input ASV table containing decimal values'),
    make_option('-d', '--divide_by', type='float',
                help='Divide ASV table values by this amount. Typical use is to "undo" the ' +
                'effect of a previous --multiply_by operation.'),
    make_option('-l', '--list_fp', type="new_filepath", help='Path to file containing a LIST' +
                'of user-supplied values, one for each ASV, with the format: <ASV><tab><value>'),
    make_option('-n', '--normalize2unity', action='store_true', help='Normalize the ASV table ' +
                'to unity. Each column will sum to 1.0'),
    make_option('-f', '--force', action='store_true',
                dest='force', help='Force overwrite of existing output file' +
                ' [default: %default]'),
]
SCRIPT_INFO['version'] = __version__



def load_asv_and_cnadjust_values(input_fp):
    """Loads the output file of the [rrna_copy_adjustments.py] script"""
    d = {}
    with open(input_fp, 'r') as input_file:
    
        for line in input_file:
        
            #split line into a list
            a = line.rstrip().split("\t")
            
            d[a[0]] = np.float64(a[1])

    return (d)


def main():
    """The main function"""
    # option_parser, opts, args = parse_command_line_parameters(**SCRIPT_INFO)
    opts = parse_command_line_parameters(**SCRIPT_INFO)[1]
    
    #validate supplied options
    if opts.input_fp == None:
        print ("*** Error *** Option '-i' is required")
        sys.exit(0)
    if opts.multiply_by != None and opts.divide_by != None:
        print ("*** Error *** Define either option '-m' or option '-d' but not both")
        sys.exit(0)
    if opts.list_fp != None and os.path.isfile(opts.list_fp) == False:
        print ("*** Error *** The file supplied with the '-l' (list_fp) option does not exist!")
        sys.exit(0)
    if opts.normalize2unity != None and (opts.multiply_by != None or opts.divide_by != None or opts.list_fp != None):
        print ("*** Error *** Option '-n' (normalize2unity) cannot be used with other math or list_fp options!")
        sys.exit(0)
    if opts.output_fp != None and opts.force == None and os.path.isfile(opts.output_fp) == True:
        print ("*** Warning *** Output file ["+opts.output_fp+"] exists! Use option '-f' to overwrite")
        sys.exit(0)
    
    #notify user of what will be done
    print ("Loading ["+opts.input_fp+"]")
    if opts.list_fp != None:
        print ("Dividing ASV table values by the values in: ["+opts.list_fp+"]")
    elif opts.output_fp != None and opts.multiply_by != None:
        print ("Multiplying ASV table values by: ["+str(opts.multiply_by)+"]")
    elif opts.output_fp != None and opts.divide_by != None:
        print ("Dividing ASV table values by: ["+str(opts.divide_by)+"]")
    elif opts.output_fp != None and opts.normalize2unity != None:
        print ("Normalizing ASV table to unity")
    elif opts.output_fp == None:
        pass
    else:
        print ("** Option Error ** Option --multiply_by OR --divide_by is required when using option --output_fp.")
        sys.exit(0)


    #load the biom table
    currBiom = load_table(os.path.abspath(opts.input_fp))

    #get the observations (asv ids)
    observ_ids = currBiom.ids(axis='observation')
    
    #get the sample ids
    sample_ids = currBiom.ids(axis='sample')
    
    #get the metadata
    observ_metadata = currBiom.metadata(axis='observation')
    
    #get the count data
    data = currBiom.matrix_data
    
    #TODO: check that ALL ASVs IN THE ASV TABLE ARE IN THE ADJUSTMENTS FILE
    # (Note that the reverse is not necessary - ie., it's okay if we're operating on a smaller/subset ASV table)
    if opts.list_fp:
        #load the copy-number adjustments file
        cn_adjustment = load_asv_and_cnadjust_values(opts.list_fp)
        
        #initialize a matrix to hold copy-number adjusted values
        cnadjusted = np.zeros(data.shape, dtype=np.float64)
        
        #perform the specified copy-number adjustment on each ASV
        for asv_id in cn_adjustment:
            if currBiom.exists(asv_id, axis="observation"):
                idx = currBiom.index(asv_id, 'observation')
                row_vals = currBiom.data(asv_id, axis='observation', dense=True)
                #we'll assume that we always want to divide by the cn_adjustment value...
                cnadjusted[idx,:] = row_vals / cn_adjustment[asv_id]
            else:
                print ("*** Error *** ASV ID '"+asv_id+"' could not be found in ["+opts.input_fp+"]")
                print ("(The adjustments file contains an ASV not found in the ASV table)")
                sys.exit(0)
        
        #put cnadjusted values into the original variable
        data = cnadjusted
    
    elif opts.normalize2unity:
        #initialize a matrix to hold normalized values
        normalized = np.zeros(data.shape, dtype=np.float64)
        
        #perform adjustments
        for sample_id in sample_ids:
            if currBiom.exists(sample_id, axis="sample"):
                idx = currBiom.index(sample_id, 'sample')
                col_vals = currBiom.data(sample_id, axis='sample', dense=True)
                npsum = np.sum(col_vals)
                if npsum > 0:
                    normalized[:,idx] = col_vals / npsum
                else:
                    print ("*** Note *** SampleID '"+sample_id+"' sums to zero!")
            else:
                print ("*** Error *** SampleID '"+sample_id+"' could not be found in ["+opts.input_fp+"]")
                sys.exit(0)
        
        #put new vals into original var
        data = normalized
        
    else:
        #perform selected mathematical operation
        if opts.multiply_by:
            #multiply ASV table by MULTIPLY_BY (note that the type remains float64)
            data = np.rint(data * opts.multiply_by)
            
        elif opts.divide_by:
            #divide ASV table by DIVIDE_BY
            data /= opts.divide_by
            
    #calculate the maximum possible subsampling depth that can be used by [multiple_rarefactions_even_depth.py]
    # WITHOUT DISCARDING ANY SAMPLES, which will be the smallest by-column sum of the integerized data matrix
    max_sampling_depth = np.amin(data.astype(int).sum(axis=0))
    
    
    #create a new biom table with the modified data, but using the original observation/sample IDs/observ_metadata
    if type(observ_metadata) == 'NoneType':
        newBiom = Table(data, observ_ids, sample_ids, table_id='ASV Table')
    else:
        newBiom = Table(data, observ_ids, sample_ids, observ_metadata, table_id='ASV Table')
    
    
    #delete any existing output file (or we'll get an error when trying to overwrite the existing biom file)
    if opts.output_fp and os.path.isfile(opts.output_fp):
        print ("Attempting to delete existing biom file: ["+opts.output_fp+"]")
        try:
            os.remove(opts.output_fp)
            print ("Existing biom file successfully removed.")
        except OSError:
            print ("Could not delete file: ["+opts.output_fp+"]")
            pass
    
    
    #save and notify
    if opts.output_fp:
        write_biom_table(newBiom, opts.output_fp)
        exit_message = "Saved biom file to: ["+opts.output_fp+"]"
        if opts.force:
            exit_message = "Biom file ["+opts.output_fp+"] overwritten with new data"
        print(exit_message)
    
    
    #notify user of max sampling depth
    print ("MAXIMUM POSSIBLE SUBSAMPLING DEPTH (w/o any samples being discarded): ["+str(max_sampling_depth)+"]")
    
    #if table is not saved, remind user that max_sampling_depth will change if saved with a different math op later
    math_op_performed = (opts.multiply_by != None or opts.divide_by != None)
    if math_op_performed and opts.output_fp == None:
        print ("(** NOTE ** This number ASSUMES YOU WILL SAVE and use an asv_table ",
               "using the SAME multiplication or division value as given here)")



if __name__ == "__main__":
    main()