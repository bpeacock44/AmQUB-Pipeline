# AmQUB Pipeline Tutorial

This pipeline has five parts, which should be run sequentially:

[Part 1: Sequence Preprocessing](#part-1-sequence-preprocessing) 
Removes PhiX reads, generates summary statistics (e.g., how many reads per sample), and optionally converts mismatched barcodes to perfect matches.

[Part 2: Demultiplexing and Trim Stat Generation](#part-2---demultiplexing-and-trim-stat-generation) 
Demultiplexes the data and estimates how many reads you’ll retain after trimming and quality filtering.

[Part 3: Taxonomic Unit (TU) Selection](#part-3)  
Processes the data to define TUs (e.g., OTUs or ASVs) and generates the initial count table.

[Part 4: Taxonomic Assignment](#part-4)
Assigns taxonomy to TUs using BLAST against the NCBI nt database and, optionally, a Qiime2 classifier.

[Part 5: Final Processing](#part-5)  
Generates finalized count tables (with taxonomy at multiple levels) and creates a detailed summary file for accessing results.

👉 Examples of how the full pipeline can be run are provided [here.](LINK) 

🪵 Note: Each script creates a log file in the output folder. Part 4 also generates a separate BLAST-specific log. These are helpful for troubleshooting errors.

&nbsp;

## Running AmQUB: Command-Line Arguments vs Parameter Files
When running AmQUB, you need to tell the program which files to use (e.g., FASTQ, mapping file). There are two ways to do this:
1. Direct Command-Line Arguments (copy/paste into terminal)
2. Parameter Files (saved on HPCC and supplied as input)

### 1. Command-Line Arguments (Direct Input) ➡️
This means pasting everything directly into the terminal.

Example:
```sh
AmQUB_part1.sh -f ~/path/to/data/XY.fastq -p ~/path/to/data/XY_map.txt
```
Here’s what’s happening:
AmQUB_part1.sh → the part of the pipeline you’re running
The letters with a dash before them are "flags" - these indicate the different options you are giving the code.
-f data/XY.fastq → tells AmQUB which FASTQ file to use
-p data/XY_map.txt → tells AmQUB which mapping file to use

### 2. Parameter Files (File Input) 🔄
Instead of typing everything into the command, you can save your settings in a comma-delimited file (CSV). Each line specifies one option.

Templates are provided in the [templates folder](LINK). You can edit them in a text editor or in Excel (but save as CSV).

Please keep in mind that your options cannot contain extra commas! So if your file path or column name etc. has a comma in it, the parameter file option will not work. 

Example params.csv:
```
Raw Fastq File,~/path/to/data/XY.fastq
Mapping File,~/path/to/data/XY_map.txt
```
Here, the text before the comma is the option (e.g., “Raw Fastq File”), and the text after is the input. Always provide full file paths.

⚠️ There are some options that are not required (optional options). You can either leave these rows out entirely or provide the value DEFAULT if you do not want to use them.

Run it like this:
```sh
AmQUB_part1.sh params.csv
```

### 🤨 Which Should You Use?
It’s up to you:
- Command-line arguments are easy to keep in a bash script for quick edits and re-runs.
- Parameter files may feel more approachable for beginners, and they serve as a permanent “paper trail” if you name and organize them clearly (e.g., with dates or project IDs).

&nbsp;

## Part 1: Sequence Preprocessing

### 📚 What does this part do?  
1. **Creates a new output folder**  
   - Named **`part1_XY_output`**, where `XY` is your FASTQ filename (without `.fq` or `.fastq`).  
     *Example: `XY.fastq` → `part1_XY_output`*  

2. **Removes likely PhiX reads**  
   - PhiX is a common sequencing control/contaminant. Removing it avoids false signals.  
   - Generates:  
     - **`XY.phiX_clean.fq`** → FASTQ file without PhiX reads  
     - **`XY.phiX_clean.alnout`** → Detailed alignment log of PhiX removal  

3. **Handles barcode mismatches (optional)**  
   - Can correct barcodes with up to *N* mismatches (see `-m` flag).  

4. **Creates demultiplexing input files**  
   - **Barcode-only file**: `XY_BC.M#.fq`  
   - **Processed reads file**: `XY.M#.fq` (if mismatches corrected), or a simulated link if not  

5. **Generates a barcode read count file**  
   - **`XY.M#.read_counts.txt`** → read counts per barcode  

6. **Generates a read statistics file**  
   - **`XY.M#.fastq_info.txt`** → contains letter frequencies, read length distribution, and error estimates (EE values)  

### ⚙️ Usage
This script accepts input in two ways: **direct command-line arguments** or a **parameter file**.

#### 1. Direct Command-Line Arguments
```sh
AmQUB_part1.sh -f data/XY.fq -p data/XY_map.txt -m 2
```
**Required Flags:**
- `-f`  Raw Fastq File  
- `-p`  Mapping File  

**Optional Flags:**
- `-m`  Mismatch Bases

#### 2. Using a Parameter File
```sh
AmQUB_part1.sh params.csv
```
Option Headers:
```
Raw Fastq File
Mapping File
Mismatch Bases
```
Example filled out:
```
Raw Fastq File,data/XY.fq
Mapping File,data/XY_map.txt
Mismatch Bases,2 #(OPTIONAL ROW)
```

### 🔍 In-Depth Description of Parameters

#### 1. **Raw Fastq File (-f) [Required]**
  - Must be **raw, uncompressed FASTQ files** (ending in `.fastq` or `.fq`).  
  - No trimming, filtering, or demultiplexing should have been done yet.  

#### 2. **Mapping File (-p) [Required]**
  - **Tab-delimited** file with at least these two columns:  
```
#SampleID   BarcodeSequence   
F001.236    CTCGACTACTGA    
F002.236    TGACCAGTAGTC    
F003.236    GCGATTAGGTCG    
```
The first two columns **must** be exactly as shown:  
  - Columns must be named exactly as shown.  
  - **#SampleID** → Sample identifier  
  - **BarcodeSequence** → Corresponding barcode sequence  
  - ⚠️ If using the mismatch option, **all original barcodes must be included** for accurate results.  

#### 3. **Mismatch Bases (-m) [Optional, default = 0]**  
  - Allows correction of reads with barcodes containing up to *N* mismatches.  
  - Range: **1–5**.  

**Example:**  
  - True barcode: `CTCGACTACTGA`  
  - Read barcode: `TTCGACTACTGA`  
  - With `-m 0` → excluded  
  - With `-m 1` → included  

⚠️ *Best practice:* Avoid using this option unless necessary, and keep the number of allowed mismatches low. 

### 📤 What to Do Next  

After **Part 1**, you may want to exclude certain samples before proceeding to **Part 2**.  

Some reasons to exclude samples:  
- **Biological filtering:** e.g., exclude males if analyzing only females  
- **Low quality/coverage:** remove samples with too few reads  

👉 **Make a copy of your original mapping file and delete unwanted rows** before moving to Part 2. I like to modify the original mapping file name with a subset ID. For example JB236_map.txt can be JB236_males_map.txt.  

&nbsp;

## Part 2: Demultiplexing and Trim Stat Generation

### 📚 What does this part do?  
1. **Creates a new output folder**  
   - Named **`part2_XY_output`**, following the same format as Part 1.  

2. **Generates a barcode reference file**  
   - **`barcodes.fa`** → contains the barcodes from your mapping file.  

3. **Demultiplexes reads using the barcode file**  
   - Reads with barcodes not found in the file are removed.  
   - Sample ID information is added to the read headers.  
   - Output file: **`XY.M#.demux.fq`**  

4. **Generates statistics for retained reads based on trim lengths**  
   - Output file: **`XY.M#.eestats.txt`** → shows how many reads remain at different trim lengths.  

#### 🤨 What are trim lengths?  
Sequencing read ends often have more errors, so trimming them improves quality. But trimming also shortens reads, which may reduce how easily they can be classified.  

The goal is to balance **quality** (fewer errors) and **length** (enough sequence for classification).  

Filtering will also remove low-quality reads in later steps, so trimming now increases the number of usable reads.  

### ⚙️ Usage
This script accepts input either by **direct command-line arguments** or via a **parameter file**.  

#### 1. Direct Command-Line Arguments
```sh
AmQUB_part2.sh -f part1_XY_output -p data/XY_males_map.txt -m 2 -r 200-301 -i 10 -s males
```
**Required flags:**  
- `-f` → Part 1 Output Folder  
- `-p` → Mapping File (can be a subset of the previous one)  

**Optional flags:**  
- `-m` Mismatch Bases
- `-s` Subset ID
- `-r` Trim Length Stats Range
- `-i` Trim Length Stats Interval

#### 2. Using a Parameter File
```sh
AmQUB_part2.sh params.csv
```
Option Headers:
```
Part 1 Output Folder
Mapping File
Mismatch Bases #(OPTIONAL ROW)
Subset ID #(OPTIONAL ROW)
Trim Length Stats Range #(OPTIONAL ROW)
Trim Length Stats Interval #(OPTIONAL ROW)
```
Example filled out:
```
Part 1 Output Folder,part1_XY_output
Mapping File,data/XY_males_map.txt
Mismatch Bases,2
Subset ID,males
Trim Length Stats Range,200-301
Trim Length Stats Interval,10
```

### 🔍 In-Depth Description of Parameters

#### 1. **Part 1 Output Folder (`-f`) [Required]**
The path to the Part 1 output folder for the flowcell being processed. (e.g. part1_AB_output)

#### 2. **Mapping File (`-p`) [Required]**
A tab-delimited file containing at least two columns:
```
#SampleID   BarcodeSequence   
F001.236    CTCGACTACTGA    
F002.236    TGACCAGTAGTC    
F003.236    GCGATTAGGTCG    
```
  - The headers **#SampleID** and **BarcodeSequence** must be exactly as shown.
  - **#SampleID** → Sample identifier  
  - **BarcodeSequence** → Barcode sequence  
  - If you subset this file by deleting rows, only the remaining samples will be processed.  

#### 3. **Subset ID (`-s`) [Optional, Default is no subset]**
Appends the subset ID to the output folder name.  

#### 4. **Mismatch Bases (`-m`) [Optional, Default = 0]**
  - Number of mismatches allowed in barcodes.  
  - Use the same value as in Part 1 (or don't set it if you didn't set it.)

#### 5. **Trim Length Stats Range (`-r`) [Optional, Default is full length of the read]**
  - Range of read lengths to generate statistics for.  
  - Example: `250-301` → stats for lengths 250 through 301.  

#### 6. **Trim Length Stats Interval (`-i`) [Optional, Default is every 25 and each of the highest 5]**
  - Interval for generating stats.  
  - Example: `5` with range `250-301` → stats for lengths 250, 255, 260 … 300.  

Note that increasing the interval and range will increase the time it takes to calculate the stats!

### 📤 What Should I Do with the Output?
After **Part 2**, check the stats output file: **`XY.M#.eestats.txt`**.  

**Example Output:**
```
16737356 reads, max len 301, avg 300.8

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
   270    8398084( 50.2%)   10343240( 61.8%)   11992976( 71.7%)
   280    7929776( 47.4%)   10008096( 59.8%)   11775144( 70.4%)
   290    7361744( 44.0%)    9599660( 57.4%)   11535620( 68.9%)
   300    6523644( 39.0%)    8902236( 53.2%)   11061544( 66.1%)
```
This file shows how many reads remain at different trim lengths under varying stringency levels (**MaxEE** values).  

- **MaxEE** (Maximum Expected Errors) works like a p-value → lower is more stringent.  
- Use the table to decide on a trim length:  
  - If merging multiple flowcells in Part 3, all must use the same trim length.  
  - Pick a length that balances read retention with high quality.  

👉 If you realize you need different stats, rerun Part 2 with adjusted parameters before moving on.  

&nbsp;

## Part 3: Taxonomic Unit (TU) Selection

### 🤨 What is a TU?  
A **Taxonomic Unit (TU)** is a general term we use to describe the “bins” or “groups” of sequences in your dataset. Each TU represents a distinct biological entity. If you see more reads in one TU than another, that suggests the corresponding entity is more abundant.  

TUs can be defined in different ways, depending on the algorithms and assumptions used. In this pipeline you can choose to work with **ASVs** or **OTUs**, which are two specific kinds of TUs. Unless describing steps unique to one method, I’ll refer to them generically as TUs.  

🦠 **ASVs vs. OTUs**  
- **OTUs (Operational Taxonomic Units):** Created by clustering sequences based on similarity (e.g., all sequences ≥97% similar to a reference or “seed” sequence are grouped into the same OTU).  
- **ASVs (Amplicon Sequence Variants):** Represent exact sequences after careful error correction. Algorithms detect and remove sequencing errors and chimeras (artificial PCR artifacts). This makes ASVs more precise and reproducible than OTUs.  
**Best practice:** Compare results from both OTUs and ASVs. If the overall taxonomic composition is similar, that builds confidence in your conclusions.  


Best practice: Many researchers compare results from both OTUs and ASVs. If the overall taxonomic composition is similar between the two, that builds confidence in your conclusions.

### 📚 What does this part do?  
1. **Creates a new output folder**  
   - Named according to your parameters.  

2. **Trims reads**  
   - Trims the demultiplexed reads from Part 2 for each included folder.  
   - Output: **XY_###bp.fq** (one per folder).  

3. **Creates a read pool fastq file**  
   - Combines trimmed reads into one file: **combined.fq**.  
   - Used later for TU assignment.  

4. **Filters trimmed reads**  
   - Filters reads using **MaxEE** (default 1.0 or specified).  
   - Output per folder: **XY_###bp.filtered.fa**.  
   - Combined output: **filtered.fa**.  

5. **Finds unique sequences in the filtered reads**  
   - Identifies unique sequences in **filtered.fa**.  
   - Removes low-abundance sequences (threshold set in parameters).  
   - Output: **uniques.fa**.  

6. **Picks TU seed sequences**  
   - Runs **UNOISE3 (ASVs)** or **UPARSE (OTUs)** on **uniques.fa**.  
   - Output: **asvs.fa** or **otus.fa** in corresponding subfolders.  

7. **Creates a TU table**  
   - Sorts reads from **combined.fq** into TU bins.  
   - Outputs table of read counts per TU and per sample: **asv_table_00.txt** or **otu_table_00.txt**.  

8. **Post-processing**  
   - Sorts TU tables.  
   - Produces a **.biom** table for use in other programs.  
   - Annotates **otus.fa** or **asvs.fa** headers with counts.  
   - Moves annotated fasta into a **blast** folder for Part 4.  

9. **Alternative optional strategies**  
   - **Strategy 1 (default):** All samples pooled first.  
   - **Strategy 2:** TU seeds generated separately for groups of samples (defined by metadata columns), then pooled and deduplicated. Helps avoid high-abundance samples skewing TU assignments.  
   - **Strategy 3:** Bin new reads into pre-existing TU seeds from a previous analysis. Useful for comparing results across projects.  
   - Outputs from these strategies are stored under `otus/STRATEGY2` or `otus/STRATEGY3`.  

**More Explanation About Alternative Strategies 2 and 3**
    - **Strategy 2** follows the same process but with a slight reordering: TUs are created separately for groups of samples instead of pooling all samples together. Usually, all samples are pooled first, but this can be problematic if your samples are very different. We developed Strategy 2 while processing soils from different states. Pooling all samples sometimes caused distinct groups to be lumped together, especially if one sample type had very high abundance of certain organisms. By processing samples separately, we reduce the risk of TU assignments being skewed by highly abundant sequences. After the TU seed sequences are picked, they are pooled and de-duplicated before creating the final TU table.
    - **Strategy 3** is simpler and useful if you want to compare TUs between two separate analyses. Sometimes TUs differ because they were processed separately. In this case, you can bin your new reads into TUs defined by a previous analysis, using the seed sequences from that earlier analysis.

### ⚙️ Usage
This script accepts input in two ways: **direct command-line arguments** or a **parameter file**.

#### 1. Direct Command-Line Arguments
```sh
AmQUB_part3.sh --in part2_AB_males_output+part2_XY_males_output --len 200 --out output_dir --al UPARSE --fr yes --fe yes --fm 0.5 --fd 2 --ug yes --us 2 --min 8 --ag yes --pid 0.98 --map merged_map.txt --col Soil_Type+Tissue --pre IDx_OTUs.fa --tblid 0.98 --un yes --rmf yes   
```

Note that these flags have two dashes before them! Usually, if a flag is one letter it only needs one dash. If it's multiple letters, it needs two. (There are exceptions.)

**Required Flags:**
- `--in` → Part 2 output folders to process (joined with `+`)  
- `--len` → Trim length  
- `--out` → Output folder name  
- `--al` → Algorithm (UPARSE for OTUs, UNOISE3 for ASVs) 

**Optional Flags:**
- `--fr` FILTERING-Generation of Removed Reads File
- `--fe` FILTERING-EE Value Appended to IDs
- `--fm` FILTERING-MaxEE Value
- `--fd` FILTERING-Discard if there are > n Ns in the read
- `--ug` UNIQUES-Generation of Sequence Binning Report
- `--us` UNIQUES-Set minuniquesize
- `--min` ALGORITHM-Set Min Size
- `--ag` ALGORITHM-Generation of Processing Report
- `--pid` UPARSE-Set Percent Identity Threshold
- `--alpha` UNOISE-Set Alpha Parameter
- `--map` STR2-Mapping File 
- `--col` STR2-Treatment Column
- `--pre` STR3-Pre-existing TUs 
- `--tblid` TABLE-Minimum fractional ID
- `--un` TABLE-Generation of Fastq of Unassigned Seqs
- `--rmf` TABLE-Generation of Read Mapping File

#### 2. Using a Parameter File
```sh
AmQUB_part3.sh params.csv
```
Option Headers:
```
Part 2 Output Folders To Process Together
Trim Length
Output Folder
UPARSE or UNOISE3
FILTERING-Generation of Removed Reads File #(OPTIONAL ROW)
FILTERING-EE Value Appended to IDs #(OPTIONAL ROW)
FILTERING-MaxEE Value #(OPTIONAL ROW)
FILTERING-Discard if there are > n Ns in the read #(OPTIONAL ROW)
UNIQUES-Generation of Sequence Binning Report #(OPTIONAL ROW)
UNIQUES-Set minuniquesize #(OPTIONAL ROW)
ALGORITHM-Set Min Size #(OPTIONAL ROW)
ALGORITHM-Generation of Processing Report #(OPTIONAL ROW)
UPARSE-Set Percent Identity Threshold #(OPTIONAL ROW)
UNOISE-Set Alpha Parameter #(OPTIONAL ROW)
STR2-Mapping File #(OPTIONAL ROW)
STR2-Treatment Column #(OPTIONAL ROW)
STR3-Pre-existing TUs #(OPTIONAL ROW)
TABLE-Minimum fractional ID #(OPTIONAL ROW)
TABLE-Generation of Fastq of Unassigned Seqs #(OPTIONAL ROW)
TABLE-Generation of Read Mapping File #(OPTIONAL ROW)
```
Example filled out (for UPARSE - note that the UNOISE3 parameter is omitted):
```
Part 2 Output Folders To Process Together,part2_AB_males_output+part2_XY_males_output
Trim Length,200
Output Folder,output_dir
UPARSE or UNOISE3,UPARSE
FILTERING-Generation of Removed Reads File,true
FILTERING-EE Value Appended to IDs,true
FILTERING-MaxEE Value,0.5
FILTERING-Discard if there are > n Ns in the read,2
UNIQUES-Generation of Sequence Binning Report,true
UNIQUES-Set minuniquesize,2
ALGORITHM-Set Min Size,8
ALGORITHM-Generation of Processing Report,true
UPARSE-Set Percent Identity Threshold,0.98
STR2-Mapping File,merged_map.txt
STR2-Treatment Column,Soil_Type+Tissue
STR3-Pre-existing TUs,IDx_OTUs.fa
TABLE-Minimum fractional ID,0.98
TABLE-Generation of Fastq of Unassigned Seqs,true
TABLE-Generation of Read Mapping File,true
```

### 🔍 In-Depth Description of Parameters

Note that most of these commands utilize USEARCH. If you want a more in-depth look at what the parameters do, please go read the documentation! It's very helpful.

#### 1. **Part 2 Output Folders To Process Together (`--in`) [Required]**
The path to each of the the Part 2 output folders you want to process with a "+" between each one. (e.g. part2_AB_male_output+part2_XY_male_output)

#### 2. **Trim Length `(--len)` [Required]**
Indicate the trim length you want to use. (See the section on Part 2 for more in-depth discussion of trim length selection.) This will be applied to all the datasets indicated in the -in flag.

#### 3. **Output Folder `(--out)` [Required]**
This can be whatever you want! Name it something that helps you remember which datasets you combined. (If I am running a lot of combinations, I might just name one "table_1_output, table_2_output" and then create a key elsewhere to remember what each table is.)

#### 4. **UPARSE or UNOISE3 `(--al)` [Required]**
Here you indicate which algorithm you want to use, which determines if you're using OTUs or ASVs. OTUs are UPARSE and ASVs are UNOISE3.

The rest of the parameters are OPTIONAL.

<ins>Parameters for the read filtering section of the pipeline:</ins>
Read more about this command [here.](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### 5. **Generation of Removed Reads File `(--fr)` [Optional, Default is to not generate it.]**
This option will generate a file that contains the reads that were removed when filtering. It will be called **`XY_###bp_discarded.fq`** where ###bp is the trim length. 

#### 6. **EE Value Appended to IDs `(--fe)` [Optional, Default is to not append them.]**
This appends the EE value (quality score) to each read ID.

#### 7. **MaxEE Value `(--fm)` [Optional, Default is 1.0.]**
This is the MaxEE filtering value. All reads that fall above this value will be discarded.

#### 8. **Discard if there are > n Ns in the read `(--fd)` [Optional, Default is no removal.]**
This option will discard any reads that contain more Ns than the number you indicate.

<ins>Parameters related to the command that finds unique reads in your filtered dataset.</ins>
Read more about this command [here.](https://drive5.com/usearch/manual/cmd_fastx_uniques.html)

#### 9. **Generation of Sequence Binning Report `(--ug)` [Optional, Default is to note generate it.]**
This file will show you which representative sequence each read was binned under. It will be called **uniques_binning_report.txt**.

#### 10. **Set minuniquesize `(--us)` [Optional, Default is 1.]**
If you set this option, any sequences that have fewer copies than the number you indicate will be discarded. (By default everything is kept because there is at least 1 copy of each sequence.)

<ins>Parameters related to TU seed sequence picking.</ins>
Read more about UNOISE3 (ASVs) [here.](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
Read more about UPARSE (OTUs) [here.](https://www.drive5.com/usearch/manual/cmd_cluster_smallmem.html)

#### 11. **Set Min Size `(--min)` [Optional, Default is 8.]**
This pecifies the minimum abundance. Input sequences with lower abundances are discarded.

#### 12. **Generation of Processing Report `(--ag)` [Optional, Default is to not generate it.]**
This report gives you more detail on the processing for each sequence. It will be called **unoise3_processing_report.txt** or **uparse_processing_report.txt**.

#### 13. **UPARSE ONLY - Set Percent Identity Threshold `(--pid)` [Optional, Default is 0.97.]**
This is threshold of similarity after which sequences will be considered to be similar enough to be of the same OTU. (So if you put 0.97, sequences that are 97% similar or more will be grouped together.)

#### 14. **UNOISE ONLY - Set Alpha Parameter `(--alpha)` [Optional, Default is 2.0.]**
The -unoise_alpha option specifies the alpha parameter (see [UNOISE2 paper](https://www.drive5.com/usearch/manual/citation.html) for definition - it's more messy but basically it's determining how stringent or lenient the algorithm is in determining what is and isn't an ASV and finding chimeras.)

<ins>Parameters for using Strategy 2.</ins> If you set these, Strategy 2 will be implemented. If you don't, it will be skipped.

#### 15. **Mapping File `(--map)` [Optional, Default is no Strategy 2.]**
This is a tab-delimited combined mapping file with all the samples from all the datasets that you are combining for this TU table. To make it, you can just merge all the mapping files from each of the folders you indicated in the "--in" flag. Make sure to remove any duplicate headers so you only have one! The #sampleID should be kept the same so your reads can still map to the correct sample. You should also have a column (or multiple) that show the groupings you want to use! For example:

```sh
#SampleID   BarcodeSequence SoilNumber
F025.236    AGACGGTAACAT    3           
F034.239    ACGGATGCAGAA    3
F073.236    GTGCTAGCTGTT    5
F037.239    TAGAGTGGTCGC    5
F014.236    TGACACACTGAA    7
F040.239    GTTGTGTTAGCC    7
F062.236    GACATGCTAATC    9
F043.239    CTATCAGACTGG    9
F039.236    CGATCCTCCAGA    10
F046.239    TGCTTATCGTCT    10
```
In this case I want to group my samples by SoilNumber. So the two samples from SoilNumber 3 will be pooled and TU seed sequences will be picked from that pool. Same with all the other Soil Numbers. Then each soil's TU seed sequences will be pooled together and "deduplicated" (identical sequences removed) before moving onto mapping the reads to the TUs.

#### 16. **Treatment Column (`--col`) [Optional, Default is no Strategy 2.]**
This should be the column name as it appears in your mapping file. (e.g. SoilNumber as above). 

<ins>Parameters for using Strategy 3.</ins> If you set this, Strategy 3 will be implemented. If you don't, it will be skipped.

#### 17. **Pre-existing TUs (`--pre`) [Optional, Default is no Strategy 3.]**
This is a fasta file of pre-existing TU seed sequences generated in another analysis that you want to map your reads to.

<ins>Parameters for making the TU table:</ins>
Read more about this command [here.](https://www.drive5.com/usearch/manual/cmd_otutab.html)

#### 18. **Minimum fractional ID (`--tblid`) [Optional, Default is 0.97 which is like 97%.]**
This is the % similarity that reads share with an TU seed sequence in order to be mapped to that TU. If you are running UPARSE, you likely want this to match the Percent Identity Threshold which is 0.97 by default.

#### 19. **Generation of Fastq of Unassigned Seqs (`--un`) [Optional, Default is to not generate one.]**
This file contains all the sequences that were not assigned to an TU. It will be called unbinned_reads.fq. 

#### 20. **Generation of Read Mapping File (`--rmf`) [Optional, Default is to not generate one.]**
This file shows which TU each read was mapped to. It will be called read_binning_report.txt.

### 📤 What Should I Do with the Output?
After running Part 3, your **TU table is ready**. The next step is to assign taxonomy (Part 4).  
- Use **BLAST** to compare TU sequences against the NCBI nt database.  
- Optionally, also run a **Qiime2 classifier** for a complementary taxonomy assignment.  

&nbsp;

## Part 4: Taxonomic Assignment 

### 📚 What does this part do?  
Part 4 assigns **taxonomy** to the TUs (OTUs or ASVs) created in Part 3. You can do this via **BLAST** searches against NCBI’s nt database, a **QIIME2 classifier**, or both. You can also choose which strategy outputs from Part 3 to annotate.  
1. **Generates a BLAST command**  
    - Creates a script (**`blast.sh`**) for running BLAST.  

2. **Runs BLAST**  
    - It generates up to 1000 matches or "hits" for each TU.
    - It checks whether 1000 hits were enough to see a decrease in bitscore.
        - If not, it does another BLAST on just those TUs to produce 30k hits. (This is rare.)
    - Checks completion and ensures **`final.blastout`** is produced for each. This and all the following BLAST-related files are found in the otus/blast or asvs/blast folder. 

3. **Filters the final.blastout file for environmental sample hits**
    - This is done by generating a list of "likely" environmental samples based on taxonomic key words. Any hit with a taxonomic ID in the list is removed.
    - Output file is called **`filtered.blastout`**

4. **Parses the BLAST output into a summary file keeping only the top bitscore hits**
    - Output is called **`top_hit_summary.txt`**

5. **Assigns LCA to each TU using the summary file**
    - All hits that matched with the top bitscore are considered. The least common taxonomic level between all of them is used as the assignment.
    - For example, if there were these 3 hits, all with the same bitscore: 
        1. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia
        2. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Klebsiella
        3. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Shigella
        - The assignment would be Enterobacteriaceae.
    - Another example: 
        1. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia
        2. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae
        3. Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; Enterobacterales
        - The assignment would be Enterobacterales. 
    - Output is called **`tax_assignments.txt`**

6. Counts how many taxa were identified down to each taxonomic level
    - Output is called **`taxa_levels.txt`** 

7. **Optionally assigns taxonomy via classifier**  
    - Uses a QIIME2 classifier to assign taxonomy in addition to BLAST.
    - Output file is found in the **`classifier_output`** folder in your asvs/otus folder. It's called **`taxonomy.tsv`**

8. **Creates a normalized version of the ASV table.**
    - This is called **`asv_table_01.norm.txt`** or **`otu_table_01.norm.txt`**

9. **Creates a summary file of taxonomic assignments** 
    - File is called **`taxonomy_summary_table.tsv`**
    - The file contains the taxonomic assignments from BLAST (taxonomy column) and the classifier (classifier_taxonomy column) if used 
    - Also includes the percent ID (blast_per_ID) and percent query coverage (blast_per_qcov) for the top BLAST hit
    - Includes confidence of classifier assignment (classifier_confidence) 
    - Shows the average abundance of each TU (avg_abun) and gives the sequence (sequence)     
    - Contains the following columns: 
        - ID (ID for the TU)
        - taxonomy (BLAST-derived taxonomic assignment)
        - blast_per_ID (the percent identification of the top BLAST hit)
        - blast_per_qcov (the percent query coverage of the top BLAST hit)
        - classifier_taxonomy (QIIME2-classifier-derived taxonomic assignment)
        - classifier_confidence (confidence level of the classifier assignment)
        - avg_abun (average abundance of the TU calculated from normalized counts)
        - sequence (the sequence for the TU)
   

### ⚙️ Usage  

Two input options: **direct command-line arguments** or a **parameter file**.  

#### 1. Direct Command-Line Arguments  

Required flags:  
- `-i` → Part 3 Output Folder To Process  
- `-t` → Number of Threads Available   
- `-e` → Email

Optional flags:  
- `-b` → Type of BLAST 
- `-v` → Expect Threshold for BLAST 
- `-s` → Skip BLAST
- `-c` → Classifier 
- `-f` → Confidence Interval 
- `--nostrategy1` → Skip Assigning Taxonomy for Strategy 1
- `--strategy2` → Assign Taxonomy for Strategy 2
- `--strategy3` → Assign taxonomy for Strategy 3  

#### 2. Using a Parameter File
```sh
AmQUB_part4.sh params.csv
```
Option Headers:
```
Part 3 Output Folder To Process 
Number of Threads Available 
Email 
Type of BLAST #(OPTIONAL ROW)
Expect Threshold for BLAST #(OPTIONAL ROW)
Skip BLAST #(OPTIONAL ROW)
Classifier #(OPTIONAL ROW)
Confidence Interval #(OPTIONAL ROW)
Skip Assigning Taxonomy for Strategy 1 #(OPTIONAL ROW)
Assign Taxonomy for Strategy 2 #(OPTIONAL ROW)
Assign Taxonomy for Strategy 3 #(OPTIONAL ROW)
```
Example filled out: 
```
Part 3 Output Folder To Process,output_dir
Number of Threads Available,256
Email,email@email.com
Type of BLAST,megablast
Expect Threshold for BLAST,0.005
Skip BLAST,true
Classifier,classifier.qza
Confidence Interval,disable
Skip Assigning Taxonomy for Strategy 1,true
Assign Taxonomy for Strategy 2,true
Assign Taxonomy for Strategy 3,true
```

### 🔍 In-Depth Description of Parameters  
 
#### 1. **Part 3 Output Folder To Process (`-i`) [Required]**  
   - This is the output directory created in part 3.
   - It must contain the TUs from Part 3 (`otus` or `asvs`).  
   - STRATEGY2 and STRATEGY3 subdirectories required if you want to assign taxonomy to them.  

#### 2. **Number of Threads Available (`-t`) [Required]**  
   - Used to speed up BLAST and classifier computations. 
   - See [01_setup_tutorial.md](LINK) for information on requesting resources for this part.

#### 3. **Email (`-e`) [Required]**  
   - Required by NCBI for BLAST queries. Must be a valid email address.  

The rest of the parameters are OPTIONAL.

#### 4. **Type of BLAST (`-b`) [Optional, default is blastn]**
   - This is the type of BLAST that will be run.  

#### 5. **E-value (`-v`) [Optional, default is 0.001]**
   - Maximum expect value for BLAST hits.  

#### 6. **Skip BLAST (`-s`) [Optional, default is no skip]** 
   - Skips the BLAST portion (though it will still process the BLAST output.)
   - Useful if you already ran BLAST but ran into an issue with the steps after.

#### 7. **Classifier (`-c`) [Optional, default is no classifier step]** 
   - Path to the classifier you want to use. See [Getting Classifiers for Qiime2](LINK) below.

#### 8. **Confidence (`-f`) [Optional, default is 0.7]** 
   - Minimum confidence for classifier assignments OR you can put "disable" to disable it entirely.

#### 9. **Strategy flags (`--nostrategy1`, `--strategy2`, `--strategy3`) [Optional, default is to run strategy 1 and not 2 or 3]** 
   - You can use these to skip assigning taxonomy to the default strategy OTUs (strategy 1) or to enable taxonomy assignment for strategy 2 or 3 respectively.

### Getting Classifiers for Qiime2
To be written. (LINK)

### 📤 What Should I Do with the Output?  
   - You can compare strategies using the **`taxonomy_summary_table.tsv`** file to determine which gives the most reliable taxonomic resolution for your dataset. 
   - If you notice any issues or you want to adjust anything you can re-run part 4.
   - You can also decide if there are TUs you don't want to keep. Make files of TUs you want to remove with one TU ID per line to use in part 5 for this purpose. You can make these for Strategy 2 and Strategy 3 separately as well.

&nbsp;

## Part 5: Final Processing 

### 📚 What does this part do?  
1. **Adds taxonomic assignments to the TU table**  
    - Will use BLAST taxonomy by default or classifier if using the "Classifier Assignments Primary" option (`-c`)
    - Table will be named **`otu_table_03_add_taxa.txt`** or **`asv_table_03_add_taxa.txt`**

2. **Generates 3 more taxonomic levels of TU tables**
    - Merge TUs at different taxonomic levels:
        - Phylum level (**`otu_table_03_add_taxa_L2.txt`** or **`asv_table_03_add_taxa_L2.txt`**)
        - Genus level (**`otu_table_03_add_taxa_L6.txt`** or **`asv_table_03_add_taxa_L6.txt`**)
        - Species level (**`otu_table_03_add_taxa_L7.txt`** or **`asv_table_03_add_taxa_L7.txt`**)

3. **Optionally, splits the TU table into 3 major domains (k__Archaea" "k__Bacteria" "k__Eukaryota)**
    - This is done with the `-u` flag and is generally only used in cases of universal assays where a universal primer was used.
    - AFTER splitting the raw count table, normalized tables will be generated as well.

4. **Generates a TU table with the sequences added**
    - Also contains taxonomic assignments.
    - Output table is called **`otu_table_04_add_seqs.txt`** or **`asv_table_04_add_seqs.txt`**

5. **Optionally, conducts the mixed top 10 analysis**
    - This step checks the consistency of the top 10 BLAST hits for each TU:
        - If all top hits belong to the same family (or higher taxonomic group) → the TU is marked as not mixed.
        - If the top hits span different families → the TU is flagged as mixed.
    - A mixed flag in the output doesn’t mean the TU is invalid, but it does indicate that its taxonomy is ambiguous. This usually happens if:
        - The sequence is too short for confident assignment.
        - The organism is underrepresented in the reference database.
        - The sequence comes from a highly conserved region shared by many taxa.

6. Optionally, generates a detailed summary file 
    - This will be called **`Detailed_Informational_otu_Table.tsv`** or **`Detailed_Informational_asv_Table.tsv`**
    - Contains the following columns: 
        - ID (ID for the TU) 
        - Sequence (the sequence for the TU)     
        - taxonomy (BLAST-derived taxonomic assignment)
        - per_ID (the percent identification of the top BLAST hit)
        - per_qcov (the percent query coverage of the top BLAST hit)        
        - c_taxonomy (QIIME2-classifier-derived taxonomic assignment)
        - confidence (confidence level of the classifier assignment)      
        - avg_abun  (average abundance of the TU calculated from normalized counts)    
        - mixed (optional - results of the mixed top 10 analysis mentioned in #5 above. Yes = mixed. No = not mixed.)
        - The rest of the normalized TU table follows (regular columns named by samples as before)
    - The first row has "nd" in all of these columns and then for each sample contains a sum of the raw counts for each sample. This can be helpful to check how many reads were assigned to each TU.

7. All tables will be generated in three formats - .txt (tab-delimited), .biom (BIOM format), and .qza (Qiime2 compatible) 

### ⚙️ Usage  

Two input options: **direct command-line arguments** or a **parameter file**.  

#### 1. Direct Command-Line Arguments  

Required flags:  
- `-i` → Part 3 Output Folder To Process 

Optional flags:  
- `-1` → Taxonomic Units to Keep File 
- `-2` → Taxonomic Units to Keep File STRATEGY2
- `-3` → Taxonomic Units to Keep File STRATEGY3
- `-u` → Universal Assay
- `-c` → Classifier Assignments Primary  
- `-d` → Detailed Informational TU Table
- `-m` → Mixed Top 10 Analysis 
- `--nostrategy1` → Skip Strategy 1 
- `--strategy2` → Process Strategy 2 
- `--strategy3` → Process Strategy 3 

#### 2. Using a Parameter File
```sh
AmQUB_part5.sh params.csv
```
Option Headers:
```
Part 3 Output Folder To Process
Taxonomic Units to Keep File #(OPTIONAL ROW)
Taxonomic Units to Keep File STRATEGY2 #(OPTIONAL ROW)
Taxonomic Units to Keep File STRATEGY3 #(OPTIONAL ROW)
Universal Assay #(OPTIONAL ROW)
Classifier Assignments Primary #(OPTIONAL ROW)
Detailed Informational TU Table #(OPTIONAL ROW)
Mixed Top 10 Analysis #(OPTIONAL ROW)
Skip Strategy 1 #(OPTIONAL ROW)
Process Strategy 2 #(OPTIONAL ROW)
Process Strategy 3 #(OPTIONAL ROW)
```
Example filled out: 
```
Part 3 Output Folder To Process,output_dir
Taxonomic Units to Keep File,to_keep.txt
Taxonomic Units to Keep File STRATEGY2,to_keep_S2.txt
Taxonomic Units to Keep File STRATEGY3,to_keep_S3.txt
Universal Assay,true
Classifier Assignments Primary,true
Detailed Informational TU Table,true
Mixed Top 10 Analysis,true
Skip Strategy 1,true
Process Strategy 2,true
Process Strategy 3,true
```

### 🔍 In-Depth Description of Parameters  
 
#### 1. **Part 3 Output Folder To Process (`-i`) [Required]**  
   - This is the output directory created in part 3 that has also been run through part 4.
   - STRATEGY2 and STRATEGY3 subdirectories required if you want to assign taxonomy to them.  

The rest of the parameters are OPTIONAL.

#### 2. **Taxonomic Units to Keep File (`-1`) [Optional, default is no removal]**  
   - This is a text file. Each line is a TU ID. Any TUs in this file will be removed from all output tables.

#### 3. **Taxonomic Units to Keep File STRATEGY2 (`-2`) [Optional, default is no removal]**  
   - This is a text file. Each line is a TU ID. Any TUs in this file will be removed from all output tables.

#### 4. **Taxonomic Units to Keep File STRATEGY3 (`-3`) [Optional, default is no removal]**
   - This is a text file. Each line is a TU ID. Any TUs in this file will be removed from all output tables.

#### 5. **Universal Assay (`-u`) [Optional, default is no splitting]**
   - If you run a universal assay, this option can be used to split your results by domain before normalizing and further processing.
   - These will be k__Archaea, k__Bacteria, and k__Eukaryota

#### 6. **Classifier Assignments Primary (`-c`) [Optional, default is to use the BLAST assignments]** 
   - When creating the TU tables, the assignments generated by the Qiime2 classifier will be used instead of the BLAST-derived assignments.

#### 7. **Detailed Informational TU Table (`-d`) [Optional, default is to not generate this table]** 
   - This enables the creation of the Detailed Informational TU Table. (See details above under [📚 What does this part do?](LINK]))

#### 8. **Mixed Top 10 Analysis (`-m`) [Optional, default is to not do this analysis]** 
   - This causes the mixed top 10 analysis to be conducted. (See details above under [📚 What does this part do?](LINK)

#### 9. **Strategy flags (`--nostrategy1`, `--strategy2`, `--strategy3`)** 
   - You can use these to skip processing the default strategy tables (strategy 1) or to enable processing of the strategy 2 or 3 tables respectively.

### 📤 What Should I Do with the Output?  
You're finished! Your otus/asvs folder (in the part3 output directory) is populated with all the tables in various formats. At this point there are optional analyses you can move forward with that are detailed in the next tutorial. You can also use Qiime2 tools to analyze further.

   


