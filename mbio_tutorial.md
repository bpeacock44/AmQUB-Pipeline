# Tutorial

This pipeline is in four parts. Each part should be run sequentially.

[Part 1:](#part-1) This first part will take as input a fastq file and a mapping file (or a batch of them). It will check for barcode mismatches and can optionally be used to convert mismatches to perfect matches up to the number of mismatches designated by you. (This is only recommended if you have samples with very few reads and you need more but if that is the case, you may want to resequence regardless!)

[Part 2:](#part-2) This part will use usearch to generate basic stats about your fastq file (which you should take a look at) and generate stats about the effect that trimming and filtering will have on your read counts. (e.g. if you trim to x length, then you will have x percent of your reads left over after the filtering step occurs.) You will need to determine a length here to use in the next part. 

[Part 3:](#part-3) This part will use usearch trim and filter your reads in preparating for ASV picking. (Note that filtered reads will be used for ASV picking, but all your reads will be used in the final ASV table.) It will then pick your ASVs and create your ASV table. 

[Part 4:](#part-4) This part is optional, as you may want to use other methods for assigning taxonomy to your ASVs. But Part 4 uses BLAST to assign taxonomy to your ASVs using the NCBI nt database. This is useful for any kind of targeted analyses where you do not have a good curated database available for assignment. (e.g. a universal assay) This pipeline can be manipulated in various ways to try to improve taxonomic assignments, given that the NCBI nt database is not curated. Rather than assigning taxonomy using the top hit, it finds all hits that have the highest bitscore and finds the least common ancestor (LCA) among them, giving preference to non-environmental samples and taxonomic groups that you may or may not choose to indicate.

There are some [examples](#example-of-overall-pipeline) of how this might be run overall at the end. 

**Note** that when you run any of these scripts, a log file will be created in you working directory that you can reference if you run into errors. Please raise an issue if this occurs! There will be a blast log created separately from your part 4 log specific to the blast run and it can be found in your final analysis output folder.

## Part 1
### USAGE
This script expects to be given at least two aguments:
- -d: a working directory, which contains a one folder for each of your fastq files named by ID
- -j: all the IDs you intend to process in a comma-delimited list (ID1,ID2,ID3,etc.)

Optional argument:
- -m: number of mismatched bases (OPTIONAL - if you want to convert barcode with given number of mismatches into perfect match barcodes)

Examples:
```sh
mbio_part1.sh -d /path/to/dir -j ID1,ID2 
mbio_part1.sh -d /path/to/dir -j ID1,ID2 -m 1
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

IMPORTANT NOTE: At this stage, the map file should contain the barcodes for ALL SAMPLES present in the fastq file or the results may contain undetected errors.

## Part 2

### USAGE
This script expects to be given at least two arguments:
- -d: a working directory, which contains the folder containing the fastq file you want to process.
- -j: a single ID. This script must be run individually on your IDs 
(This is in contrast to part 1, which was run just once for all.)

Optional arguments:
- -m: the number of mismatches you want to use. This needs to match the files you generated in part 1.
- -s: this is a comma-delimited list of trim lengths you want to view stats for. These will be generated in addition to all cutoffs that are 11 bases or fewer below the max length of your reads.
- -o: a subset ID - if you want to run further analyses on a subset of the samples in your data, you can create a mapping file in the same format as the original with the lines of unwanted samples removed. This file will be named ID_map.subsetID.txt (e.g. JB141_map.Nickels01.txt) and be placed in the same ID folder as the other files are in.

Examples:
```sh
mbio_part2.sh -d /path/to/dir -j ID1 
mbio_part2.sh -d /path/to/dir -j ID1 -s 137,148
mbio_part2.sh -d /path/to/dir -j ID1 -o sub1 
mbio_part2.sh -d /path/to/dir -j ID1 -m 1
mbio_part2.sh -d /path/to/dir -j ID1 -o sub1 -m 1
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
This script expects to be given at least 4 arguments:
- -d: a working directory, which contains one folder for each of your fastq files named by ID
- -j: the folders created in the last part that you intend to process in a comma-delimited list (ID1_subset1_output, ID2_output, ID3_subset2_output, etc.)
- -l: the length you want to trim your reads to. Note ALL files will be trimmed to this length.
- -o: the name of your output directory

Optional arguments:
- -m: number of mismatches, if using (again, this should have been specified from part1)
- -n: change the minsize of the unoise3 algorithm (default is 8)

Examples:
```sh
mbio_part3.sh -d /path/to/dir -j "ID1_output,ID2_output" -l 150 -o test1_out
mbio_part3.sh -d /path/to/dir -j "ID1_output,ID2_output" -l 150 -o test2_out -m 1
mbio_part3.sh -d /path/to/dir -j "ID1_output,ID2_output" -l 150 -o test3_out -n 8
```

### INPUT ###
This script follows part 2, which must be completed first. You will look at your trim stats, determine what length you want to trim to, and run this code to finish the analysis.

So, as an example, your working directory might now include:
- Folder JB141_Nickels01_output and directory JB143_output, both containing output of part2.

When this code is run, a new directory named as you indicated will be created for analysis output. 

## Part 4

This part can be run a few different ways, depending on how you want to run BLAST. 

1) Run within the singularity as an interactive session on a cluster (request the resources you will need for BLAST)
2) Run locally on your computer within the singularity
3) If you have a very large ASV file, you may want to split it up and run on separate resources. In this case, you will have a slightly more complicated pipeline. See the "Part 4 Split" section below. For both 1 and 2, the following applies.

### USAGE 
This script expects to be given at least 5 arguments:
- -d: a working directory, which contains one folder for each of your fastq files named by ID
- -o: the name of your output directory
- -b: the path to your blast script file
- -e: email of the user for NCBI purposes

Optional arguments:
- -t This is a filter file - it will have four tab-delimited columns. 
Include any taxonomic groups that you want to preferentially keep or reject as well as their taxonomic ID. 
ALL taxonomies included under these taxonomic IDs will be treated accordingly so check NCBI
and make sure. If you are doing a universal assay, do not include the -t flag and DO include the -u flag.
- -u: universal assay - causes final ASV tables to be split into taxonomic groups prior to normalizing
- -s: skip the blast - skips the blast portion - useful for troubleshooting or re-running taxonomy assignment steps etc. Note that if -s is enabled, -b is not required.
- -j: this flag creates a specialized excel summary output that Dr. Borneman specifically requested. Runtime will increase, as it requires an analysis examining the top 10 blast hits for each ASV.

Examples:
```sh
mbio_part4.sh -d /path/to/dir -o test1_out -b /path/to/blast.sh -e email@email.com -t filterfile.txt 
mbio_part4.sh -d /path/to/dir -o test3_out -e email@email.com -s
mbio_part4.sh -d /path/to/dir -o test4_out -b /path/to/blast.sh -e email@email.com -u -j
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
blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 
```
If running a local blast, you could use the same format but do not include the SBATCH and module load lines.

For argument -t, you can choose to include a filter file that will essentially indicate taxonomic groups you want to give preference for and reject outright in the taxonomic assignment process. For example, if I was running the analysis on bacterial ITS amplicon data taken a plant sample, then I might want to give preference to bacterial taxa and reject any plant taxa.

The columns will be as follows, with a header: (Name, ID, Rank, and Action)
- Name - the taxonomic name on NCBI
- ID - the taxonomic ID on NCBI
- Rank - the rank of the taxonomic group (Note that this doesn't have to match NCBI, as seen here with bacteria, which is listed as a superkingdom on NCBI. It is primarily for your information.)
- Action - whether you want to give preference to (keep) or reject this group.

```
Name    ID    Rank    Action
Bacteria    2    k    Keep
Viridiplantae    33090    k    Reject
```
Note that rank is lowercase.

When assigning taxonomy, decisions will be made based on bitscore. Highest bitscore results among non-rejected taxonomic groups will contribute to the taxonomic assignment. If the highest bitscores are in your "keep" group or the the highest bitscore is shared between groups, then your taxonomic assignment will be from the "keep" group. However, if the bitscore is higher in environmental samples of the "keep" group OR, secondarily, if it is higher in unspecified groups (taxonomies not listed in your filter file) then those will be used instead.

## Part 4 Split

If you want to split up your ASV file, you will need to run the blast portion on it's own outside of the container. You will start with the script "mbio_part4_SPLIT_blast.sh", which is a truncated version of mbio_part4_blast.sh.

It has no optional arguments:
- -d: a working directory, which contains one folder for each of your fastq files named by ID
- -o: the name of your output directory
- -b: the path to your blast script file (should include ALL jobs, each beginning with a shebang as below)
- -n: number of different jobs you want to run
- -r: the type of blast run you want to do (local or slurm)

### INPUT 
As before, this script follows part 3, which must be completed first. The output directory will have already been generated in part 3.

For argument -b, you are going to want to make a blast script but it should include ALL the scripts you want to run. The code will automatically split this file by shebang, and also split your ASV file and merge the results back together into a single blast file.
```
# ########## SLURM BLAST FILE EXAMPLE (-b), if you had 2 jobs to run ########## 
#!/bin/bash 
#SBATCH -p epyc
#SBATCH -c 256
#SBATCH --mem=300G 
# any other parameters or modules needed
module load ncbi-blast/2.6.0+
module load db-ncbi

#<>#<>#<>#<>#<>
# YOU MUST SET THESE:
#<>#<>#<>#<>#<>
DATABASE_PATH=$NCBI_DB/nt
NUMTHREADS=256

#<>#<>#<>#<>#<>
# GENERALLY DON'T CHANGE THESE:
#<>#<>#<>#<>#<>
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
TASK=blastn
INFASTA=$1
MAXTSEQS=$2  
EVAL=0.001
blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 

#!/bin/bash 
#SBATCH -c 256
#SBATCH --mem=300G 
# any other parameters or modules needed
module load ncbi-blast/2.6.0+
module load db-ncbi

#<>#<>#<>#<>#<>
# YOU MUST SET THESE:
#<>#<>#<>#<>#<>
DATABASE_PATH=$NCBI_DB/nt
NUMTHREADS=256

#<>#<>#<>#<>#<>
# GENERALLY DON'T CHANGE THESE:
#<>#<>#<>#<>#<>
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
TASK=blastn
INFASTA=$1
MAXTSEQS=$2  
EVAL=0.001
blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 
```

After you run mbio_part4_SPLIT_blast.sh, you will need to run mbio_part4_blast.sh as usual - just make sure to run with the -s flag so it doesn't run BLAST all over again!

## Example of Overall Pipeline:

```sh
# start an interactive node with some power.
srun -c 128 --mem 300gb --pty bash -l
module load singularity

#A directory called sing_to_bind contains usearch and targeted_microbiome_via_blast, which is cloned from this repository UNTIL ADDED TO SINGULARITY
#WDIR contains sing_to_bind and the sif file, as well as your folders with data.

## Enter the singularity using this command. 
WDIR=/path/to/WDIR
singularity shell --bind ${WDIR}/sing_to_bind:/bind/ ${WDIR}/1.2.2.sif --cleanenv --no-home 

WDIR=/path/to/WDIR
SDIR=${WDIR}/sing_to_bind

##### pre-processing, barcode mismatches, etc. (Can run on any number of files.) 
${SDIR}/targeted_microbiome_via_blast/mbio_part1.sh -d ${WDIR} -j "ID1"

## INTEGRATE new part2 - doesn't seem to be in the github cloud yet. Improve instructions for this as well
##### demultiplexing, collecting stats. (Must be run on individual files AFTER part 1a has been run.)
${SDIR}/targeted_microbiome_via_blast/mbio_part2.sh -d ${WDIR} -j "ID1" -o small
${SDIR}/targeted_microbiome_via_blast/mbio_part2.sh -d ${WDIR} -j "ID2" -o big

##### trimming, filtering, combining, ASV generation, ASV table generation, etc. in preparation for taxonomic classification
${SDIR}/targeted_microbiome_via_blast/mbio_part3.sh -d ${WDIR} -j "ID1_output,ID2_output" -l 300 -o PN1_final_results

# must add targeted_microbiome_via_blast to the PATH UNTIL ADDED TO SINGULARITY
export PATH="${SDIR}/targeted_microbiome_via_blast/helper_functions:$PATH"
${SDIR}/targeted_microbiome_via_blast/mbio_part4_blast.sh -d ${WDIR} -o PN1_final_results -e beth.b.peacock@gmail.com -b blast.sh -r slurm

```

## Example of Overall Pipeline with Split Blast:

```sh
# start an interactive node with some power.
srun -c 128 --mem 300gb --pty bash -l
module load singularity

#A directory called sing_to_bind contains usearch and targeted_microbiome_via_blast, which is cloned from this repository.
#WDIR contains sing_to_bind and the sif file, as well as your folders with data.

## Enter the singularity using this command. 
WDIR=/path/to/WDIR
singularity shell --bind ${WDIR}/sing_to_bind:/bind/ ${WDIR}/1.2.2.sif --cleanenv --no-home 

WDIR=/path/to/WDIR
SDIR=${WDIR}/sing_to_bind

##### pre-processing, barcode mismatches, etc. (Can run on any number of files.) 
${SDIR}/targeted_microbiome_via_blast/mbio_part1.sh -d ${WDIR} -j "ID1"

## INTEGRATE new part2 - doesn't seem to be in the github cloud yet. Improve instructions for this as well
##### demultiplexing, collecting stats. (Must be run on individual files AFTER part 1a has been run.)
${SDIR}/targeted_microbiome_via_blast/mbio_part2.sh -d ${WDIR} -j "ID1" -o small
${SDIR}/targeted_microbiome_via_blast/mbio_part2.sh -d ${WDIR} -j "ID2" -o big

##### trimming, filtering, combining, ASV generation, ASV table generation, etc. in preparation for taxonomic classification
${SDIR}/targeted_microbiome_via_blast/mbio_part3.sh -d ${WDIR} -j "ID1_output,ID2_output" -l 300 -o PN1_final_results

#EXIT singularity AND node bc to run blast you will be submitting two jobs.
exit
exit 

WDIR=/path/to/WDIR
SDIR=${WDIR}/sing_to_bind

module load db-ncbi # this must be available locally
module load seqkit # this must be available locally

# must add targeted_microbiome_via_blast to the PATH
export PATH="${SDIR}/targeted_microbiome_via_blast/helper_functions:$PATH"
${SDIR}/targeted_microbiome_via_blast/mbio_just_blast.sh -d ${WDIR} -o PN1_final_results -b blast_to_split.sh -n 2 -r slurm

srun -c 128 --mem 300gb --pty bash -l
module load singularity
WDIR=/path/to/WDIR
singularity shell --bind ${WDIR}/sing_to_bind:/bind/ ${WDIR}/1.2.2.sif --cleanenv --no-home 
WDIR=/path/to/WDIR
SDIR=${WDIR}/sing_to_bind

##### classification and final file generation
${SDIR}/targeted_microbiome_via_blast/mbio_part4_blast.sh -d ${WDIR} -o PN1_final_results -e beth.b.peacock@gmail.com -s

```
