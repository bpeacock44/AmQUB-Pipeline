# Tutorial

This pipeline is in four parts. Each part should be run sequentially.

[Part 1:](#part-1) This first part will take as input a fastq file and a mapping file (or a batch of them). It will check for barcode mismatches and can optionally be used to convert mismatches to perfect matches up to the number of mismatches designated by you. (This is only recommended if you have samples with very few reads and you need more but if that is the case, you may want to resequence regardless!)

[Part 2:](#part-2) This part will use usearch to generate basic stats about your fastq file (which you should take a look at) and generate stats about the effect that trimming and filtering will have on your read counts. (e.g. if you trim to x length, then you will have x percent of your reads left over after the filtering step occurs.) You will need to determine a length here to use in the next part. 

[Part 3:](#part-3) This part will use usearch trim and filter your reads in preparating for ASV picking. (Note that filtered reads will be used for ASV picking, but all your reads will be used in the final ASV table.) It will then pick your ASVs and create your ASV table. 

[Part 4:](#part-4) This part is optional, as you may want to use other methods for assigning taxonomy to your ASVs. But Part 4 uses BLAST to assign taxonomy to your ASVs using the NCBI nt database. This is useful for any kind of targeted analyses where you do not have a good curated database available for assignment. (e.g. bacterial ITS) Rather than assigning taxonomy using the top hit, it finds all hits that have the highest bitscore and finds the least common ancestor (LCA) among them.

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
Each folder needs to contain a fastq file named by ID followed by "\_raw.fq" and an appropriately formatted mapping file followed by "\_map.txt"
For example, the directory indicated contains a folder called JB141 and it contains the files "JB141_raw.fq" and "JB141_map.txt."

Your mapping file should be formatted in the same way that the Qiime program uses:

At least three tab-delimited columns:
1) SampleID (with a # before as in #SampleID) - these are the IDs you associate with each sample
2) BarcodeSequence - the barcodes for each sample

Other columns can include any characteristics you want to use later on to run differential analysis or correlation, etc.
```
SampleID	BarcodeSequence	SampleType	PlatePosition	Library	TubeLabel	Contents	DateTaken
B001.110	CTCGACTACTGA	SAMPLE	A1	JB110	1	Psyllid	1-6	2/28/19
B002.110	TGACCAGTAGTC	SAMPLE	A2	JB110	2	Psyllid	7-12	2/28/19
B003.110	GCGATTAGGTCG	IGNORE	A3	JB110	3	Psyllid	13-18	2/28/19
PCR_CONTROL	ACATGGCCTAAT	CONTROL	A4	JB110	NA	NA	NA
```

IMPORTANT NOTE: At this stage, the map file should contain the barcodes for ALL SAMPLES present in the fastq file or the results may contain undetected errors.

If you are confused about the mapping file, there are some more notes [here](#more-mapping-file-details).

## Part 2

### USAGE
This script expects to be given at least two arguments:
- -d: a working directory, which contains the folder containing the fastq file you want to process.
- -j: a single ID. This script must be run individually on your IDs 
(This is in contrast to part 1, which was run just once for all.)

Optional arguments:
- -m: the number of mismatches you want to use. This needs to match the files you generated in part 1.
- -s: this is a comma-delimited list of trim lengths you want to view stats for. These will be generated in addition to all cutoffs that are 11 bases or fewer below the max length of your reads.
- -o: a subset ID - if you want to pick ASVs from and further analyze only a subset of the samples in your data, you can create a mapping file in the same format as the original with the lines of unwanted samples removed. This file will be named ID_map.subset_name.txt (e.g. JB141_map.florida.txt) and be placed in the same ID folder as the other files are in.

Examples:
```sh
mbio_part2.sh -d /path/to/dir -j ID1 
mbio_part2.sh -d /path/to/dir -j ID1 -s 137,148
mbio_part2.sh -d /path/to/dir -j ID1 -o sub1 
mbio_part2.sh -d /path/to/dir -j ID1 -m 1
mbio_part2.sh -d /path/to/dir -j ID1 -o sub1 -m 1
```
### INPUT
This script can only be run once the original fastq file (e.g. JB141_raw.fq) has been run through part 1. 

Each folder needs to contain the fastq files resulting from part 1, which are named by ID followed by .M#.fq and \_BC.M#.fq, as well as a mapping file (either the original or the subset.)

So, as an example, your working directory might now include:
- Folder JB141 (containing JB141.M0.fq, JB141_BC.M0.fq, and JB141_map.Nickels01.txt)
- JB141_map.txt should also be present in folder if subset map isn't used.

When this code is run, a new directory will be created for your output named either with the name of your subset (e.g. JB141_Nickels01_output), or with your regular ID if no unique map was provided (e.g. JB141_output) 

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
mbio_part3.sh -d /path/to/dir -j "ID1_output,ID2_output" -l 150 -o test3_out -n 7
```

### INPUT ###
This script follows part 2, which must be completed first. You will look at your trim stats, determine what length you want to trim to, and run this code to finish the analysis.

When this code is run, a new directory named as indicated by -o will be created for analysis output. 

## Part 4

This part involves using BLAST and it can be run in two different ways:

1) Within the continer. (If you are on a cluster, you will want to request resources before you begin the container!)
3) Initially outside of the container (if you have a very large ASV file and you want to split it up and run each part in parallel to speed up the process) and then finishing inside the container. In this case, you will have a slightly more complicated pipeline. See the "Part 4 Split" section below for details on this.

### USAGE 
This script expects to be given at least 5 arguments:
- -d: a working directory, which contains one folder for each of your fastq files named by ID
- -o: the name of your output directory
- -e: email of the user for NCBI purposes

Optional arguments:
- -b: path to your blast script file (required if you don't use the -s option)
- -s: skip the blast - skips the blast and uses blast output already there (required if you don't use the -b option)

- -u: universal assay - causes final ASV tables to be split into 3 domains
- -j: create a detailed summary file with extra detail for assessing taxonomic assignments. Runtime will increase slightly, as it generates an analysis examining the top 10 blast hits for each ASV.
- -c: path to a Qiime2 classifier file you want to use for generating alternative taxonomic assignments

Examples:
```sh
mbio_part4.sh -d /path/to/dir -o test1_out -b /path/to/blast.sh -e email@email.com
mbio_part4.sh -d /path/to/dir -o test2_out -e email@email.com -s
mbio_part4.sh -d /path/to/dir -o test3_out -b /path/to/blast.sh -e email@email.com -u -j -c unite_ver10_99_04.04.2024-Q2-2024.5.qza
```

### INPUT 
This script follows part 3, which must be completed first. The output directory will have already been generated in part 3.

The NCBI nt database needs to be bound to the singularity container. This is described in [00_singularity_instructions.md](https://github.com/bpeacock44/targeted_microbiome_via_blast/blob/main/00_singularity_instructions.md) and shown in the [examples](#example-of-overall-pipeline) below. If you bind it to /database as shown in the example below, then the path indicated here should work.

For argument -b, you are going to want to make a blast script. You will need to modify the NUMTHREADS below to match the number of threads you have available (whether on your local computer requested in your interactive session). This needs to be executable.
```
#!/bin/bash
DATABASE_PATH=/database/nt 
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

## Part 4 Split (note: mbio_part4_blast_only.sh is not written yet)

NOTE: You will need [ncbi blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [seqkit](https://bioinf.shenwei.me/seqkit/) installed on your system to run this, since you will not be using the container.

You will also want to have the scripts ["mbio_part4_blast_only.sh"](https://github.com/bpeacock44/targeted_microbiome_via_blast/blob/main/pipeline_scripts/mbio_part4_blast_only.sh), ["reblast_check.pl"](https://github.com/bpeacock44/targeted_microbiome_via_blast/blob/main/helper_functions/reblast_check.pl) and ["blast_iterator.sh"](https://github.com/bpeacock44/targeted_microbiome_via_blast/blob/main/helper_functions/blast_iterator.sh) in your path and executable. 

If you want to split up your ASV file, you will need to run the blast portion on it's own outside of the container. You will start with the script "mbio_part4_blast_only.sh", which is a truncated version of mbio_part4.sh.

It has no optional arguments:
- -d: a working directory, which contains one folder for each of your fastq files named by ID
- -o: the name of your output directory
- -b: the path to your blast script file (should include ALL jobs, each beginning with a shebang as below)
- -n: number of different jobs you want to run

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

After you run mbio_part4_blast_only.sh, you will need to run mbio_part4.sh as usual - just make sure to run with the -s flag so it doesn't run BLAST all over again!

## Example of Overall Pipeline:

```sh
### WDIR contains folders with your data and your blast script.
WDIR=/path/to/WDIR
cd ${WDIR}

# Start an interactive session with enough power to run BLAST if you are in a cluster environment.
srun -c 128 --mem 300gb --pty bash -l
module load singularity

# You are binding two directories to your container here - the one containing usearch and the one containing your ncbi nt database.
# This command will start the singularity container.
singularity shell --bind /path/to/usearch:/bind/ --bind /path/to/ncbi_database:/database/ /path/to/container.sif --cleanenv --no-home 

# Define your WDIR path. It should be the directory you start in.
WDIR=$(pwd)

# part 1
mbio_part1.sh -d ${WDIR} -j "ID1"

# part 2
mbio_part2.sh -d ${WDIR} -j "ID1"
mbio_part2.sh -d ${WDIR} -j "ID2"

# part 3
mbio_part3.sh -d ${WDIR} -j "ID1_output,ID2_output" -l 300 -o PN1_final_results

# part 4
mbio_part4.sh -d ${WDIR} -o PN1_final_results -e email@email.com -b ${WDIR}/blast.sh

```

## Example of Overall Pipeline with Split Blast:

```sh
### WDIR contains folders with your data and your blast script.
WDIR=/path/to/WDIR
cd ${WDIR}

# start an interactive session with some power. You won't run BLAST yet so no need to go crazy here.
srun -c 8 --pty bash -l
module load singularity

# you are binding only one directories to your container here - the one containing usearch.
# this command will start the singularity container.
singularity shell --bind /path/to/usearch:/bind/ /path/to/container.sif --cleanenv --no-home 

# define your WDIR path. It should be the directory you start in.
WDIR=$(pwd)

# part 1
mbio_part1.sh -d ${WDIR} -j "ID1"

# part 2
mbio_part2.sh -d ${WDIR} -j "ID1"
mbio_part2.sh -d ${WDIR} -j "ID2"

# part 3
mbio_part3.sh -d ${WDIR} -j "ID1_output,ID2_output" -l 300 -o PN1_final_results

# exit singularity 
exit

# exit your interactive session if you want more resources for BLAST
exit 

# make sure WDIR is still defined
WDIR=/path/to/WDIR

# save the scripts "mbio_part4_blast_only.sh", "reblast_check.pl", and "blast_iterator.sh" to WDIR and make them executable.

# add WDIR to path so these scripts can be run.
export PATH="${WDIR}:${PATH}"

# Load modules required if you are in a cluster environment. 
module load seqkit 

# part 4 (just blast)
mbio_part4_blast_only.sh -d ${WDIR} -o PN1_final_results -b blast_to_split.sh -n 2 

# start another interactive session with some power
srun -c 60 --pty bash -l

# start singularity as before
module load singularity
WDIR=/path/to/WDIR
cd ${WDIR}
singularity shell --bind /path/to/usearch:/bind/ /path/to/container.sif --cleanenv --no-home 

# define your WDIR path. It should be the directory you start in.
WDIR=$(pwd)

# part 4 (the rest, but skip blast with -s flag!)
mbio_part4.sh -d ${WDIR} -o PN1_final_results -e email@email.com -s
```

## More Mapping File Details
When you demultiplex (part 1), mapping files must only contain the samples (and barcodes) from that specific fastq file. (And you want ALL the samples to be there as well so it demultiplexes properly!) So each data file you get should have it's own mapping file.

The only time you alter a mapping file during the pipeline is if you want to make an ASV table from only a subset of the samples within a fastq file, which you can indicate at part 2 using the -o option. (For example, if you were sequencing the microbiomes of insects and your advisor asked if he could add a couple of mouse gut samples into the library as well, you probably don't want to include them - especially when you are picking ASVs! So in step 2 you would create a new mapping file without those mouse gut sample rows.)

There is no reason to combine mapping files until you are doing analyses at the very end when you will be referencing metadata and your ASV table contains samples from multiple fastq files. At this point, it makes sense to merge all the rows from all the mapping files that represent the samples in the ASV table you will be analyzing.

