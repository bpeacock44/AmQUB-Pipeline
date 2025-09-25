# Pipeline Tutorial

This pipeline is in five parts. Each part should be run sequentially.

[Part 1:](#part-1---sequence-preprocessing) Sequence Preprocessing (Removes PhiX reads, generates stats about the data and how many reads per sample you have, and optionally can convert mismatched barcodes to perfect match barcodes.)

[Part 2:](#part-2---demultiplexing-and-trim-stat-generation) Demultiplexing and Trim Stat Generation (Demultiplexes the data and generates stats about how many reads you will retain based on trim length and quality filtering.)

[Part 3:](#part-3) Taxonomic Unit (TU) Selection (Processes the data in order to pick taxonomic units and generates initial count table.)

[Part 4:](#part-4) Taxonomic Assignment (Uses BLAST of the NCBI nt database and optionally a classifier from Qiime2 to assign taxonomy to TUs.)

[Part 5:](#part-5) Final Processing (Generates different versions/levels fo the count tables with taxonomy and creates a Detailed Summary File for accessing results.)

There are some [examples](#example-of-overall-pipeline) of how this might be run overall at the end. 

**Note** that when you run any of these scripts, a log file will be created in the output folder that you can reference if you run into errors. Please raise an issue if this occurs! There will be a blast log created separately from your part 4 log specific to the blast run. 

&nbsp;

## Part 1 - Sequence Preprocessing

### 📚 What is this part of the pipeline doing?
1. **Creates a new folder** where all the output for this part will be stored.
    - This will be called **“part1_XY_output”**, where **XY** is the name of your input fastq file with `.fq` or `.fastq` removed (e.g., `XY.fastq` generates `part1_XY_output`).
2. **Removes reads that are likely phiX**
    - This step will generate:
      - **`XY.phiX_clean.fq`** – The new data file without phiX reads.
      - **`XY.phiX_clean.alnout`** – A detailed output file for the alignment process used to find phiX reads.
3. **Finds and converts mismatches in the read barcodes (if desired)**
    - See below for details on what this is and how it works.
4. **Creates the two files required to demultiplex**
    - **Barcode-only file**: Formatted like the original fastq file, named **`XY_BC.M#.fq`**.
    - **Processed reads file**: If mismatched barcodes are converted, it will create **`XY.M#.fq`**; otherwise, it creates a simulated link.
5. **Generates a barcode read count file**
    - This will be called **`XY.M#.read_counts.txt`**.
6. **Generates a read statistics file**
    - This file will be called **`XY.M#.fastq_info.txt`** and includes stats such as letter frequencies, min/med/max read lengths, and mean/med/max EE values.

### ⚙️ Usage
This script accepts input in two ways: **direct command-line arguments** or a **parameter file**.

#### 1. Direct Command-Line Arguments
```sh
AmQUB_part1.sh -f data/AB.fq -p data/AB_map.txt -m 2
AmQUB_part1.sh -f data/XY.fq -p data/XY_map.txt
```
**Required Flags:**
> **-f**  Path to fastq file for your flowcell  
> **-p**  Path to mapping file  

**Optional Flags:**
> **-m**  Number of allowed mismatches (1-5). Default is **0**.

#### 2. Using a Parameter File
```sh
AmQUB_part1.sh params.csv
```
The parameter file (e.g., `params.csv`) should be comma-delimited, without headers or extra spaces. The first column should contain the following row names:
```
Raw Fastq File
Mapping File
Mismatch Bases
```
Example with Two Flowcells:
```
Raw Fastq File,data/AB.fq,data/XY.fq
Mapping File,data/AB_map.txt,data/XY_map.txt
Mismatch Bases,2,DEFAULT #(OPTIONAL ROW)
```
- `DEFAULT` in **Mismatch Bases** means it will use the default (0 mismatches).
- Easily created in Excel: enter data in columns, then save as a CSV file.

### 🔍 In-Depth Description of Parameters

#### 1. **Raw Fastq File (-f) [Required]**
These are **raw** data files (before demultiplexing, filtering, or trimming). They must be **uncompressed** and end in `.fastq` or `.fq`.

#### 2. **Mapping File (-p) [Required]**
**Must be tab-delimited** and contain the following columns:
```
#SampleID   BarcodeSequence   
F001.236    CTCGACTACTGA    
F002.236    TGACCAGTAGTC    
F003.236    GCGATTAGGTCG    
```
The first two columns **must** be exactly as shown:  
- **#SampleID**: Sample ID associated with the barcode.
- **BarcodeSequence**: The barcode sequence.
**If using the mismatch function, you MUST include all original barcodes.** If any are missing, the function's results will not be trustworthy.

#### 3. **Mismatch Bases (-m) [Optional, Default = 0]**
This optional function **converts mismatched barcodes** into perfect match barcodes. You can specify **1-5 mismatches**.

Example:
- If `F001.236` has barcode **CTCGACTACTGA**, a read with **TTCGACTACTGA** would normally be excluded.
- If **1 mismatch is allowed**, it would be included.

⚠ **We recommend avoiding this function unless necessary and keeping mismatches minimal.**

### 📤 What Should I Do with the Output?
After running **Part 1**, you will decide which samples to include for **Part 2** by modifying the map file:

- **Example 1:** Exclude non-relevant samples (e.g., only analyze female nematodes and remove soil samples).
- **Example 2:** Remove samples with too few reads for analysis.

🔹 **Create a copy of the original map file and remove unwanted samples before proceeding to Part 2.**

&nbsp;

## Part 2 - Demultiplexing and Trim Stat Generation

### 📚 What is this part of the pipeline doing?
1. **Creates a new folder** where all the output for this part will be stored.
    - This will be called **“part2_XY_output”**, following the same naming format as Part 1.
2. **Generates a barcode reference file**
    - This file, **`barcodes.fa`**, contains the barcodes from your new mapping file.
3. **Demultiplexes the data based on the barcodes file**
    - Reads with barcodes **not found** in the barcodes file are removed.
    - Sample ID information is added to the header of each read.
    - The resulting fastq file is named **`XY.M#.demux.fq`**.
4. **Generates statistics for retained reads based on trim lengths**
    - Outputs a file called **`XY.M#.eestats.txt`**, detailing how many reads remain at various trim lengths.

### ⚙️ Usage
This script accepts input in two ways: **direct command-line arguments** or a **parameter file**.

#### 1. Direct Command-Line Arguments
```sh
mbio_part2.sh -f part1_AB_output -p data/AB_subset_map.txt -m 2 -r 200-301 -i 10
mbio_part2.sh -f part1_XY_output -p data/XY_map.txt
```

**Required Flags:**
> **-f**  Path to the output folder from Part 1  
> **-p**  Path to the mapping file  

**Optional Flags:**
> **-m**  Number of allowed mismatches (1-5). Default is **0**.  
> **-r**  Trim length stats range. Default is full read length.  
> **-i**  Trim length stats interval. Default is **25**.

#### 2. Using a Parameter File
```sh
mbio_part2.sh params.csv
```

The parameter file (e.g., `params.csv`) should be comma-delimited, without headers or extra spaces. The first column should contain these row names:
```
Part 1 Output Folder
Mapping File
Mismatch Bases
Trim Length Stats Range
Trim Length Stats Interval
```

Example with Two Flowcells:
```
Part 1 Output Folder,part1_AB_output,part1_XY_output
Mapping File,data/AB_subset_map.txt,data/XY_map.txt
Mismatch Bases,2,DEFAULT #(OPTIONAL ROW)
Trim Length Stats Range,200-301,DEFAULT #(OPTIONAL ROW)
Trim Length Stats Interval,10,DEFAULT #(OPTIONAL ROW)
```
- `DEFAULT` will use the script's default values.
- Each column represents a different flowcell.

### 🔍 In-Depth Description of Parameters

#### 1. **Part 1 Output Folder (-f) [Required]**
The path to the Part 1 output folder for the flowcell being processed. The folder name should follow the format **part1_XY_output**, where **XY** is the flowcell ID.

#### 2. **Mapping File (-p) [Required]**
A tab-delimited file containing at least two columns:
```
#SampleID   BarcodeSequence   
F001.236    CTCGACTACTGA    
F002.236    TGACCAGTAGTC    
F003.236    GCGATTAGGTCG    
```
- The headers **#SampleID** and **BarcodeSequence** must be exactly as shown.
- The first column contains sample IDs.
- The second column contains barcode sequences.
- If you subset this file to exclude unwanted samples, only those will be processed in Part 2.

#### 3. **Mismatch Bases (-m) [Optional, Default = 0]**
The number of barcode mismatches allowed. If you used this option in Part 1, use the same value here.

#### 4. **Trim Length Stats Range (-r) [Optional]**
The range of read trim lengths for which statistics should be generated. Example:
- **`250-301`** calculates stats for reads trimmed to lengths **250 through 301**.

#### 5. **Trim Length Stats Interval (-i) [Optional]**
Determines the interval for calculating trim length stats. Example:
- **`5`** calculates stats every **5 bases** within the specified range.
- If `250-301` is the range, then stats will be generated for lengths **250, 255, 260, …, 300**.

### 📤 What Should I Do with the Output?
After running **Part 2**, check the stats output file: **`XY.M#.eestats.txt`**.

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
- This table shows how many reads remain at different trim lengths based on quality filtering.
- **MaxEE** (Maximum Expected Errors) functions like a p-value: lower values are more stringent.

🔹 **Decide on a Trim Length:**
- If you plan to merge multiple flowcells in **Part 3**, all should be trimmed to the same length.
- Compare stats across flowcells and choose a trim length that retains enough reads while ensuring quality.

You can rerun **Part 2** with different parameters if needed before moving on.

# TUTORIAL TO BE CONTINUED.


