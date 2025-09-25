# Getting Set Up
NOTE: All information related to accessing services for HPCC at UCR here have been copied from their website at hpcc.ucr.edu. I will update this if I become aware of changes, but if you run into any issues with my instructions here or the links stop working, go to hpcc.ucr.edu directly to get up-to-date information. (Also, feel free to message me so I can make necessary changes.)

## Using HPCC at UCR
First you need to get a user account. Please email user account requests to support@hpcc.ucr.edu. Include the full name, NetID and email address of both users and PI. Users need to be members of the PI’s group. Preferentially, user account requests should come from the corresponding PI directly. If the request comes from a new user then the PI needs to be CC’ed in the email exchange.

Once you have an account, you will access it in one of two ways. First, you can use the terminal directly. This is the most flexible method, but you will need to learn some basic code as you will not be able to point/click on anything. If you plan to do a lot of work on HPCC, I highly recommend taking the time to learn to use the terminal. More detailed instructions are given [here](https://hpcc.ucr.edu/manuals/access/login/#a-ssh-login-from-terminal) but here are the basics:

1. Open your ssh terminal of choice. (For Linux and Apple, you can use native terminal apps but there are others with more features. On Windows, MobaXterm is a recommended.) Type the following ssh login command from a terminal application, where <username> needs to be replaced by the actual account name of a user. The <> characters indicate a placeholder and need to be removed. Next, press enter to execute the ssh command.

```sh
ssh <username>@cluster.hpcc.ucr.edu
```

2. Type your password and hit enter. Note, when typing the password the cursor will not move and nothing is printed to the screen. 

3. Follow the Duo multifactor authenication instructions printed to the screen.

## How will you access and upload your files?
The files on HPCC are accessible directly via the terminal - you can view and edit them directly using the command line. If you'd prefer, you can use OnDemand to upload/download/view the files directly in your browser. The HPCC OnDemand instance is located here: https://ondemand.hpcc.ucr.edu/. Log in with your cluster login details (UCR NetID and Password) and verify your login with Duo’s two-factor authentication. This will take you directly to your personal directory on HPCC.







## Setting Up "Working Directory"
You need to create a folder that you will be working in. (Your "working directory.") This is where you will put your input data and mapping files (in a folder called "data") and it is where your output files will be generated.



## Prep – Starting Singularity

### Steps
1. Make sure singularity is loaded as a module on the cluster (command below) or installed on your computer. Note that installing singularity on Mac or Windows will require you to use a virtual Linux machine because it is only compatible with Linux. There are detailed instructions on how to do this here (link forthcoming).
```sh
module load singularity
```
2. Create a folder of the programs/tools required. In the following examples, I have one called “mbio_pipeline_files.” And inside I have:
    - The singularity file (AmQUB.sif)
    - A copy of USEARCH 64-bit
    - A QIIME2 classifier (instructions link forthcoming)
3. Set this folder equal to the variable “programs” using the full path as so: 
```sh
programs="/path/to/mbio_pipeline_files"
```
4. You can check if step 3 was done correctly by listing the files in the folder using this command:
```sh
ls $programs
```
5. You should also have an updated local NCBI nt database (instructions) or access to one on your cluster.
6. Set your database equal to the variable “NCBI_DB” using the full path like so:
```sh
NCBI_DB=”/srv/projects/db/ncbi/preformatted/20240807”
```
7. Move to your working directory – which is where you want to generate output.
8. Run the following command to start the AmQUB singularity container. Note that there are a few lines, and they all need to be run at the same time! 
```sh
singularity shell 
--bind ${programs}:/home/programs/ \
--bind ${NCBI_DB}:/database/ \
--bind $(pwd):/home/analysis/ \
${programs}/AmQUB.sif'
```
9. Detailed explanation of what the above command means: 
- Remember that when you start a singularity container you are essentially starting a “minicomputer”, and these commands help set up what that computer has access to. You are using the flag “--bind” to set folders in your minicomputer:
    - Your programs folder will be the folder “/home/programs” in the container
    - Your NCBI database will be in the path “/database” in the container
    - Your working directory (which is stored as $(pwd) if you are in it currently) will be the folder “/home/analysis” in the container
    - Backslashes (\\) mean that the command continues on the next line.
    - The end of the command should be the singularity container file you want to start.
10. Once you run the singularity shell command, it will open the singularity and you will see “Singularity>” in your terminal. This is the prompt for you to start running commands!
```sh
Singularity>                             
```
&nbsp;