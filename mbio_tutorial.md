# Tutorial

This pipeline is in four parts. Each part should be run sequentially.

Part 1: This first part will take as input a fastq file and a mapping file (or a batch of them). It will check for barcode mismatches and can optionally be used to convert mismatches to perfect matches up to the number of mismatches designated by you. (This is only recommended if you have samples with very few reads and you need more but if that is the case, you may want to resequence regardless!) 

Part 2: This part will use usearch to generate basic stats about your fastq file (which you should take a look at) and generate stats about the effect that trimming and filtering will have on your read counts. (e.g. if you trim to x length, then you will have x percent of your reads left over after the filtering step occurs.) You will need to determine a length here to use in the next part.

Part 3: This part will use usearch trim and filter your reads in preparating for ASV picking. (Note that filtered reads will be used for ASV picking, but all your reads will be used in the final ASV table.) It will then pick your ASVs and create your ASV table.

Part 4: This part is optional, as you may want to use other methods for assigning taxonomy to your ASVs. But Part 4 uses BLAST to assign taxonomy to your ASVs using the NCBI nt database. This is useful for any kind of targeted analyses where you do not have a good curated database available for assignment. (e.g. a universal assay) This pipeline can be manipulated in various ways to try to improve taxonomic assignments, given that the NCBI nt database is not curated. Rather than assigning taxonomy using the top hit, it finds all hits that have the highest bitscore and finds the least common ancestor (LCA) among them, giving preference to non-environmental samples and taxonomic groups that you may or may not choose to indicate.

## Part 1
### USAGE
This script expects to be given at least two aguments:
- d: a working directory, which contains a one folder for each of your fastq files named by ID
- j: all the IDs you intend to process in a comma-delimited list (ID1,ID2,ID3,etc.)

Optional argument:
- m: number of mismatched bases (OPTIONAL - if you want to convert barcode with given number of mismatches into perfect match barcodes)

Examples:
```sh
mbio_part1.sh -d /path/to/dir -j JB141,JB143 
mbio_part1.sh -d /path/to/dir -j JB141,JB143 -m 1
```

### INPUT
Each folder needs to contain a fastq file named by ID followed by "\_L1P1.fq" and an appropriately formatted mapping file followed by "\_map.txt"
For example, the directory indicated contains a folder called JB141 and it contains the files "JB141_L1P1.fq" and "JB141_map.txt."

Your mapping file should be formatted in the same way that the Qiime program uses:

At least three tab-delimited columns:
1) SampleID (with a # before as in #SampleID) - these are the IDs you associate with each sample
2) BarcodeSequence - the barcodes for each sample
3) SampleType - you can use this to filter later on - e.g. removing controls before data analysis, removing samples that aren't relevant, etc.

Other columns can include any characteristics you want to use later on to run differential analysis or correlation, etc.
```
SampleID	BarcodeSequence	SampleType	PlatePosition	Library	TubeLabel	Contents	DateTaken
B001.110	CTCGACTACTGA	SAMPLE	A1	JB110	1	Psyllid 1-6	2/28/19
B002.110	TGACCAGTAGTC	SAMPLE	A2	JB110	2	Psyllid 7-12	2/28/19
B003.110	GCGATTAGGTCG	IGNORE	A3	JB110	3	Psyllid 13-18	2/28/19
PCR_CONTROL	ACATGGCCTAAT	CONTROL	A4	JB110	NA	NA	NA
```

## Part 2

### USAGE
This script expects to be given at least two aguments:
- d: a working directory, which contains the folder containing the fastq file you want to process.
- j: a single ID. This script must be run individually on your IDs 
(This is in contrast to part 1, which was run just once for all.)

Optional arguments:
- m: the number of mismatches you want to use. This needs to match the files you generated in part 1.
- o: a subset ID - if you want to run further analyses on a subset of the samples in your data, you can create a mapping file in the same format as the original with the lines of unwanted samples removed. This file will be named ID_map.subsetID.txt (e.g. JB141_map.Nickels01.txt) and be placed in the same ID folder as the other files are in.

Examples:
```sh
mbio_part2.sh -d /path/to/dir -j JB141 
mbio_part2.sh -d /path/to/dir -j JB141 -o Nickels01 
mbio_part2.sh -d /path/to/dir -j JB141 -m 1
mbio_part2.sh -d /path/to/dir -j JB141 -o Nickels01 -m 1
```
### INPUT
This script can only be run once the original fastq file (e.g. JB141_L1P1.fq) has been run through part 1, which is used to find and replace mismatched barcodes with perfect match barcodes. 

Each folder needs to contain the fastq files resulting from 1a, which are named by ID followed by \_A1P1.M#.fq and \_A1P2.M#.fq, as well as a mapping file (either the original or a subset.)

So, as an example, your working directory might now include:
- Folder JB141 (containing JB141_A1P1.M0.fq, JB141_A1P2.M0.fq, and JB141_map.Nickels01.txt)
- JB141_map.txt should also be present in folder if subset map isn't used.

When this code is run, a new directory will be created for your output named either with the unique identifier for your subset, if given (e.g. JB141_Nickels01_output), or it will be named after your regular ID if no unique map was provided (e.g. JB141_output) 

## Part 3

### USAGE
This script expects to be given at least 4 aguments:
- d: a working directory, which contains one folder for each of your fastq files named by ID
- j: the folders created in the last part that you intend to process in a comma-delimited list (ID1_subset1,ID2,ID3_subset2,etc.)
- l: the length you want to trim your reads to. Note ALL files will be trimmed to this length.
- o: the name of your output directory

Optional arguments:
- m: number of mismatches, if using (again, this should have been specified from part1)

Examples:
```sh
mbio_part3.sh -d ${MDIR} -j "JB141_Nickels01_output,JB143_output" -l 150 -o test1_out
mbio_part3.sh -d ${MDIR} -j "JB141_Nickels01_output,JB143_output" -l 150 -o test2_out -m 1
```

### INPUT ###
This script follows part 2, which must be completed first. You will look at your trim stats, determine what length you want to trim to, and run this code to finish the analysis.

So, as an example, your working directory might now include:
- Folder JB141_Nickels01_output and directory JB143_output, both containing output of part2.

When this code is run, a new directory named as you indicated will be created for analysis output. 

## Part 4

### USAGE 
This script expects to be given at least 4 aguments:
- d: a working directory, which contains one folder for each of your fastq files named by ID
- o: the name of your output directory
- b: the path to your blast script file
- r: the type of blast run you want to do (local or slurm)
- e: email of the user for NCBI purposes

Optional arguments:
- t: filter file (see below for details.) If you are doing a universal assay, do not include the -t flag and DO include the -u flag.
- m: number of mismatches, if using (again, this should have been specified from part1)
- u: just for universal assay - causes final ASV tables to be split into 3 taxonomic domains prior to normalizing
- s: skip the blast - skips the blast portion - useful for troubleshooting or re-running taxonomy assignment steps etc.

Examples:
```sh
mbio_part4.sh -d ${MDIR} -o test1_out -b ${MDIR}/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt 
mbio_part4.sh -d ${MDIR} -o test2_out -b ${MDIR}/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt -m 1 
mbio_part4.sh -d ${MDIR} -o test3_out -b ${MDIR}/blast.sh -e email@email.com -r local -s
mbio_part4.sh -d ${MDIR} -o test4_out -b ${MDIR}/blast.sh -e email@email.com -r local -m 1 -u 
```

### INPUT 
This script follows part 3, which must be completed first. The output directory will have already been generated in part 3.

For argument -b, you are going to want to make a blast script based on how you want to run blast locally or in a cluster environment. Here is an example of a cluster-based blast script:
```
#!/bin/bash 
#SBATCH -p i128
#SBATCH -c 128
# any other parameters or modules needed
module load blast-plus

#<>#<>#<>#<>#<>
# YOU MUST SET THESE:
#<>#<>#<>#<>#<>
DATABASE_PATH=/sw/dbs/blast_db_download/nt
NUMTHREADS=128

#<>#<>#<>#<>#<>
# GENERALLY DON'T CHANGE THESE:
#<>#<>#<>#<>#<>
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
TASK=blastn
INFASTA=$1
MAXTSEQS=$2  
EVAL=0.001
# blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 
```
If running a local blast, you could use the same format but do not include the SBATCH and module load lines.

For argument -t, you can choose to include a filter file that will essentially indicate taxonomic groups you want to give preference for and reject outright in the taxonomic assignment process. For example, if I was running the analysis on bacterial ITS amplicon data taken a plant sample, then I might want to give preference to bacterial taxa and reject any plant taxa.

The format will be as follows 
Header (Name, ID, Rank, and Action)
Name - the taxonomic name on NCBI
ID - the taxonomic ID on NCBI
Rank - the rank of the taxonomic group (Note that this doesn't have to match NCBI, as seen here with bacteria, which is listed as a superkingdom on NCBI. It is primarily for your information.)
Action - whether you want to give preference to (keep) or reject this group.

```
Name    ID    Rank    Action
Bacteria    2    k    Keep
Viridiplantae    33090    k    Reject
```
Note that rank is lowercase.

When assigning taxonomy, decisions will be made based on bitscore. Highest bitscore results among non-rejected taxonomic groups will contribute to the taxonomic assignment. If the highest bitscores are in your "keep" group or the the highest bitscore is shared between groups, then your taxonomic assignment will be from the "keep" group. However, if the bitscore is higher in environmental samples of the "keep" group OR, secondarily, if it is higher in unspecified groups (taxonomies not listed in your filter file) then those will be used instead.