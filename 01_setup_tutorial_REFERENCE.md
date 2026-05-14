# AmQUB Quick Setup Guide 🛠️  

## 1. Text Editor 📝
- Use any text editor (TextEdit, Notepad, Sublime, VS Code).  
- Save bash scripts as `.sh` (e.g., `run_code.sh`).  
- Add comments with `#`.  

Example:  
```sh
# Run AmQUB part 1
AmQUB_part1.sh -f data/JB236.fastq -p data/JB236_map.txt
```
## 2. HPCC Access 🔵🟡
### Log In
Option A (Terminal, recommended):
```sh
ssh <NetID>@cluster.hpcc.ucr.edu
```
Option B (Browser):
https://ondemand.hpcc.ucr.edu

###  Navigation Basics
```sh
cd folder_name    # move into folder
ls                # list files
pwd               # show current path
mkdir new_folder  # make folder
```
⚠️ Always run jobs in ~/bigdata, not home directory.
### Tmux (keep jobs running)
```sh
tmux new -s session_name
tmux attach -t session_name
tmux ls
# Detach: ctrl+b then d
```
## 3. Working Directory 🗂️
```sh
mkdir ~/bigdata/250926_test_run
mkdir ~/bigdata/250926_test_run/data
```
### Get Data
Download from core (replace credentials):
```sh
wget -r -e robots=off --no-parent \
--http-user='USERNAME' --http-password='PASSWORD' http://cluster.hpcc.ucr.edu/~genomics/USERNAME/2030/
```
Move files:
```sh
mv cluster.hpcc.ucr.edu/~genomics/USERNAME/2030/* .
rm -r cluster.hpcc.ucr.edu
```
Decompress:
```sh
pigz -d -c ~/shared/illumina_data/JB241_FC2030/Undetermined_S0_R1_001.fastq.gz > ~/bigdata/250926_test_run/data/JB241.fastq
```
### Mapping File Example (JB236_map.txt)
```
#SampleID   BarcodeSequence
F001.236    CTCGACTACTGA
F002.236    TGACCAGTAGTC
```
### File naming:
FASTQ → JB236.fastq

Map → JB236_map.txt

## 4. Start AmQUB Container 📦
```sh
module load singularity
programs="~/shared/mbio_pipeline_files/"
NCBI_DB="/srv/projects/db/ncbi/preformatted/20250707/"

cd ~/bigdata/250926_test_run

singularity shell \
--bind ${programs}:/home/programs/ \
--bind ${NCBI_DB}:/database/ \
--bind $(pwd):/home/analysis/ \
${programs}/AmQUB.sif
```
Inside container prompt:
```
Singularity>
```
🎉 You’re ready to run AmQUB commands!