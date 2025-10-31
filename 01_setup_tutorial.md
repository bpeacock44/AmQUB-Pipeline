# Getting Set Up 🛠️
In order to run AmQUB, you will need to get set up in a few ways: 
1. 📝 Find or download a text editor on your computer.
2. 💻 Get your account working and familiarize yourself with how to access files and use the terminal.
3. 🗂️ Set up your working directory and data.
4. 📦 Set up and start the AmQUB singularity container. 

&nbsp;

## 1. Text Editors 📝
In this pipeline you will be running some code. I've tried to minimize the amount that you need to run, but at the very least you need to copy and paste commands into the terminal. Whenever you do this, it's imperitive you keep track of the code you use at each step. Why?
1. You can keep track of what you've already done - this is vital if you want to troubleshoot. (Especially if you ask for help.)
2. You create a template to use for future work. (No need to re-invent the wheel.)

Every computer (Windows, Mac, Linux) comes with a basic text editor of some kind. (e.g. Mac has TextEdit). You can use these. There are nicer ones too that add features that can be helpful. For example, some text editors can color the words/letters of the code to make it easier to read. [(For example.)](https://mlfk3cv5yvnx.i.optimole.com/cb:HA53.300ea/w:auto/h:auto/q:mauto/f:best/https://www.ninjaone.com/wp-content/uploads/2024/01/what-is-bash-scripting.jpg)

I use Sublime Text, which has a lot of features but is also very basic/lightweight. You can try it for free as long as you'd like, though you will be regularly prompted to get a subscription if you keep using it. (Though there is no time limit.) It's available on Windows/Mac/Linux.  

So either figure out what your computer already has or try downloading one - either way, get one and then start copying in code that you use as you go through this tutorial. 

When running AmQUB, all code will be "Bash" which is a coding language. When you write bash code into a file and save it, the suffix should be ".sh" such as "250928_code.sh". Then when you open it in a text editor that supports text coloring, it'll automatically know how to color the text. 

One more note - feel free to add notes about the code you add to this file. When you add notes (not code) to a bash file like the one you are making, you put a "#" in front of each line of code. This indicates that it is a comment, not code. Here is an example:

```sh
# This is some code I am learning for the AmQUB pipeline.
AmQUB_part1.sh -f data/JB236.fastq -p data/JB236_map.txt
# The -f and the -p are flags. -f indicates my input file and -p indicates my mapping file.
```

Also note in the above code box - on the far right there is a "copy" function so you can just copy the code easily into your own bash file!

&nbsp;

## 2. Using HPCC at UCR 🔵🟡
NOTE: All information related to accessing HPCC at UCR here have been copied from [hpcc.ucr.edu](hpcc.ucr.edu). If you run into any issues with my instructions here or the links stop working, go there directly to get up-to-date information. (Also, feel free to message me so I can make necessary changes.)

### 2.1 Getting Your Account Running and Logging In 🖥️
To request an account:
1. Email support@hpcc.ucr.edu.
2. In your request, include:
    - Your full name
    - Your UCR NetID
    - Your email address
    - The PI’s name (and ideally CC them)
3. Requests should come from the PI. If you send it yourself, CC your PI to confirm.

Once you have an account, you will access it in one of two ways. 

<ins>Option A: Terminal Login (recommended for long-term use)</ins>
1. Open a terminal program. (Mac/Linux: Use the built-in Terminal app. Windows: Install MobaXterm.)
2. Type this command, replacing <username> with your UCR NetID:
```sh
ssh <username>@cluster.hpcc.ucr.edu
```
3. 🔑 Enter your password (nothing will show as you type—this is normal).
4. Follow the on-screen Duo two-factor authentication prompts.
5. More detailed instructions are given [here](https://hpcc.ucr.edu/manuals/access/login/#a-ssh-login-from-terminal).

<ins>Option B: Web Browser Login (OnDemand)</ins>
1. Easier for beginners. Lets you upload/download/view files in a browser.
2. Go to: [https://ondemand.hpcc.ucr.edu/](https://ondemand.hpcc.ucr.edu/)
3. Log in with your NetID and password, then approve Duo.
4. At the top menu: Files → Home Directory.
5. From here you can view, upload, and download files.

Note that you cannot run code directly from OnDemand! You can open a terminal from there, and there are some programs you can use to run code (e.g. Jupyter Notebooks, RStudio). You can also open a remote desktop to use GUI applications. The user-input parts of the AmQUB pipeline require you to copy and paste some basic code into the terminal.

### 2.2 How will you access and upload your files? 📂
The files on HPCC are accessible directly via the terminal - you can view and edit them directly using the command line. If you'd prefer, you can use OnDemand to upload/download/view the files directly in your browser. At the top of the screen you can click the menu "Files" and then "Home Directory". This will take you to your files on the cluster, which you can view, edit, download, and upload etc. via the buttons along the top.  

### 2.3 Basic Usage via the Terminal 🔳
When working in the terminal, you need to use linux commands to do anything. HPCC has a great overview [here](https://hpcc.ucr.edu/manuals/linux_basics/) and I highly recommend you familiarize yourself with these basic commands, as it will make your life much easier! I will go over a few basics here which you will definitely need to use to run AmQUB.

<ins>Understanding Where You Are (Nodes and Directories)</ins> 📍

When you log in, you’ll see something like:
```sh
[bpeacock@bluejay] ~$
```
- bpeacock → your username
- bluejay → the “head node” you’re on (these are named after birds). Nodes are computers so you're just in a "head" or "base" computer.
- ~ → your home directory

⚠️ WARNING. Do not run heavy jobs on head nodes! Use them only for browsing and light commands. Heavy analyses require requesting compute resources (covered later).

<ins>Moving Around</ins> 🔄

- cd folder_name → go into a folder
- cd ~/bigdata → go to your bigdata folder (where analyses should be run)
- ls → list the contents of the current folder
- pwd → show your current location

⚠️ WARNING. Always do analyses inside bigdata, not your home directory (home is too small and will fill up quickly).

<ins>Making/Naming Folders</ins> 📁

To make folders:
```sh
mkdir PN130_analysis
```
Avoid spaces! Use underscores instead. This makes it easier to run commands on files. Commands interpret spaces differently than we do.

<ins>TMUX</ins> 🪟

Lastly, when you log into the cluster, if you are inactive for a while, it tends to just log you out automatically. This is very frustrating - especially if you start running a command! When it logs you out, it ends anything you were running. 

There is a great workaround for this - it's called tmux. Tmux can be a little intimidating but you have to be able to use it to run this pipeline. Essentially, it creates a "session" of your computer that stays logged in - so even if you leave, it's still running in the background! I always work in a tmux session now just so I can quickly go back to what I was doing before I got off. If you request resources, those continue running too.

Start tmux:
```sh
tmux
```
Start a named session:
```sh
tmux new -s fungi_run
```
Attach to a session you left:
```sh
tmux attach -t fungi_run
```
List sessions:
```sh
tmux ls
```
Detaching requires you to use a tmux shortcut. When you are in tmux, any tmux-related commands you run require you to first type in **ctrl+b** and then whatever shortcut you want to use. For example, detach is **ctrl+b+d**. 

Tmux is very powerful - I like to create a session and run multiple windows in the same session so I can run code in one window and then also view/move files in another window. You can do this by typing **ctrl+b+c** (for create a new window) and then switch between them using **ctrl+b+n** (for next window). You can also split windows and more! The HPCC website has a much more thorough introduction [here](https://hpcc.ucr.edu/manuals/hpc_cluster/terminalide/).

&nbsp;

## 3. Setting Up "Working Directory" 🗂️
You need to create a folder that you will be working in. (This is your "working directory.") This is where you will put your input data and mapping files (in a folder called "data") and it is where your output files will be generated.

So go to bigdata (or a sub-folder if you want) and make your new folder. Name it as you please.
```sh
mkdir ~/bigdata/250926_test_run
mkdir ~/bigdata/250926_test_run/data # optionally, I like to also make a "data" folder for all input files.
```
To run AmQUB, you will need to supply two primary things - your data (an UNCOMPRESSED fastq) and mapping files.

### 3.1 Data 📊
When you get an email from the core telling you that your data is ready, you should immediately download it to the ~/shared/illumina_data/ folder **soon**. The email should provide you with a username and password to use. 

>Dear XYZ,  
>
>The sequence data from your recent sample submission has been demultiplexed and uploaded. The data is available at  (http://cluster.hpcc.ucr.edu/~genomics/USERNAME). 
> 
>Your login credentials are (replace "USERNAME" above with provided): 
> 
>Username: XYZ
>Password: password
>
>wget -r -e robots=off --no-parent --http-user='USERNAME' --http-password='PASSWORD' http://cluster.hpcc.ucr.edu/~genomics/USERNAME

(Note it says they demultiplexed the data but in our case that isn't always the case. Demultiplexing is when you split up the data by sample using the barcodes associated with them so you know which reads came from which samples, etc. More on that to come.)

You can actually plug the address they give you right into your browser (e.g. https://cluster.hpcc.ucr.edu/~genomics/jborneman/) and after giving the password, you can look through what is on there. The data is stored in a folder named after the flowcell number, which is given as the subject of the email from the core like this: "FC#2164 (JB236)". Inside you have the data as well as some quality control files. In Reports/html/ you can view how many reads were generated for each barcode.

⚠️ WARNING. The core deletes old data routinely - sometimes after just a month - so <ins>please copy your data safely to our shared folder **as soon as possible**.</ins> To do this, go onto the cluster and navigate to the ~/shared/illumina_data folder and make a folder there for your data. Name it with the JB identifier and your flowcell ID. For example:
```sh
mkdir ~/shared/illumina_data/JB241_FC2030 
cd ~/shared/illumina_data/JB241_FC2030 
```
Then download using the command given in the email with your username and password added:
```sh
wget -r -e robots=off --no-parent --http-user='XYZ' --http-password='password' http://cluster.hpcc.ucr.edu/~genomics/XYZ/2030/
```
This is sort of annoying, but when you download the data it keeps the same path so in your current folder you now have "cluster.hpcc.ucr.edu/~genomics/XYZ/2030/" instead of just the folder you wanted. Just move the contents of the flowcell folder to your current folder and then you can safely remove the cluster.hpcc.ucr.edu folder like so:
```sh
mv cluster.hpcc.ucr.edu/~genomics/XYZ/2030/* . 
rm -r cluster.hpcc.ucr.edu
```
The mv command moves or renames files. The asterisk (\*) stands for "everything" and the period stands for the folder you are currently in. So this command is saying: "move all the files in "cluster.hpcc.ucr.edu/~genomics/XYZ/2030/" to where I am currently."

The rm command deletes the folder designated and the -r flag means the command works recusrively or on everything INSIDE the folder as well. It will ask you to individually confirm deletion of each file. You can skip that by running rm -rf instead (f stands for FORCE) but this command is very dangerous - if you run it on a folder by accident you will not have a way to recover what you deleted - there is no "trash" in linux - rm deletes permanently. **So be careful!**

Now you will want to create an decompressed copy in your working directory. Look in to your newly created folder in the illumina data file (JB241_FC2030 in our example). 
```sh
ls JB241_FC2030 # putting "ls" before a path shows you what is in THAT folder, not your current one. "ls" by itself shows what is in your current directory.
'index.html?C=D;O=A'  'index.html?C=M;O=A'  'index.html?C=N;O=A'  'index.html?C=S;O=A'   JB241_S1_R1_001.fastq.gz   Stats
'index.html?C=D;O=D'  'index.html?C=M;O=D'  'index.html?C=N;O=D'  'index.html?C=S;O=D'   Reports                    Undetermined_S0_R1_001.fastq.gz
```
Lots of files there! If you used the standard microbiome barcodes in the Borneman Lab, you are mostly concerned with the one called Undetermined_S0_R1_001.fastq.gz. This is what your data is saved as! The data is stored as "Undetermined" because you haven't demultiplexed yet. So let's decompress it and send it to your working directory:
```sh
pigz -d -c ~/shared/illumina_data/JB241_FC2030/Undetermined_S0_R1_001.fastq.gz > ~/bigdata/250926_test_run/data/JB241.fastq
# pigz is a great tool for decompressing gzip files. 
# The -d flag means decompress and the -c flag tells it not to delete the compressed verion after decompressing. 
# The arrow tells it to send the decompressed output to a new file, which is designated after.
# Note that I am saving it as the new shortened identifier.
```

### 3.2 Mapping Files 🗺️
Mapping files should be tab-delimited text files. (Tab delimited just means that each "column" is separated by a tab.) Here is an example mapping file, where you have one row per sample and the associated barcode. This mapping file would be used to optionally find mismatched barcodes in part 1 and demultiplex in part 2.
```
#SampleID   BarcodeSequence   
F001.236    CTCGACTACTGA    
F002.236    TGACCAGTAGTC    
F003.236    GCGATTAGGTCG 
```
The first column (#SampleID) should be in every version of every mapping file you use with this data - that sample ID will become associated with each read in your data. The AmQUB tutorial will walk you through which kinds of mapping files you need at each step. AmQUB will use them to demultiplex your data. You can also add more columns as desired with other metadata. As you progress through the pipeline, you will add mapping files for different steps. 

You can create the mapping file on your desktop in a text editor, or you can create them in Excel: File → Save As → Tab Delimited Text (.txt). Then you can upload via OnDemand. (Or make it directly in a text editor in the terminal if you want to be fancy!)

### 3.3 Naming Conventions 🏷️
You need to name these files a particular way so AmQUB can find them. 
1. Name the fastq file with some identifier of your choice and ".fastq" after.
    - JB236.fastq 
2. Name the associated map with the **same** identifier with "\_map.txt" appended to it. 
    - JB236_map.txt

This identifier will be used to generate output folders in the first parts of AmQUB so I recommend you don't make them super long or complicated.

&nbsp;

## 4. Starting Singularity for AmQUB 📦
When you use programs on HPCC, you sometimes need to "load" them as modules. This is like loading software. Sometimes programs conflict with each other so that is why by default they aren't all loaded. Singularity is a module.
```sh
module load singularity
```
You need to have a folder of the programs/tools required. In our lab shared folder, I have one called “mbio_pipeline_files" where I keep the AmQUB-related files that you can use. Inside should be:
    - The singularity file (AmQUB.sif)
    - A copy of USEARCH 64-bit
    - A QIIME2 classifier if you want to use one later on. (Details later on in the tutorial.)

Set this folder equal to the variable “programs” using the full path as so: 
```sh
programs="~/shared/mbio_pipeline_files/"
```
Setting variables is helpful when using linux - it's a way to store long paths as a single word to use later. You can check if this was done correctly by listing the files in the folder using this command:
```sh
ls $programs
```
You also need access to the NCBi nt database. Check what the latest one is as follows:
```sh
ls /srv/projects/db/ncbi/preformatted/
# this will list all the versions of the blast database on the cluster. 
```
Set the path to the most recent database equal to the variable “NCBI_DB” using the full path like so:
```sh
NCBI_DB=”/srv/projects/db/ncbi/preformatted/20250707/"
```

Now move to your working directory – which is where you want to generate output.
```sh
cd ~/bigdata/250926_test_run
```
Run the following command to start the AmQUB singularity container. Note that there are a few lines, and they all need to be run at the same time! 
```sh
singularity shell 
--bind ${programs}:/home/programs/ \
--bind ${NCBI_DB}:/database/ \
--bind $(pwd):/home/analysis/ \
${programs}/AmQUB.sif'
```
Detailed explanation of what the above command means: 
When you start a singularity container you are essentially starting a “minicomputer”, and these commands help set up what that computer has access to. You are using the flag “--bind” to set folders in your minicomputer:
    - Your programs folder will be the folder “/home/programs” in the container
    - Your NCBI database will be in the path “/database” in the container
    - Your working directory (which is stored as $(pwd) if you are in it currently) will be the folder “/home/analysis” in the container
    - Backslashes (\\) mean that the command continues on the next line.
    - The end of the command should be the singularity container file you want to start.

Once you run the singularity shell command, it will open the singularity and you will see “Singularity>” in your terminal. This is the prompt for you to start running the AmQUB commands!
```sh
Singularity>                             
```
🎉 That means you’re inside the AmQUB container and ready to run pipeline commands! 
