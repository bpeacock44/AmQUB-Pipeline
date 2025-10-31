# Two notes here! 

# First, you can vastly simplify your code if you use variables. 
# You set variables (a group of letters such as "VAR") equal to other text (such as a path like "~/bigdata/folder1"). 
# Then any time you put ${VAR} it can stand in for the other text. Below I use the variable WDIR to stand in for my path.
# When you add variables into a path, it's a good idea to enclose the entire path in quotes. See below.

# Second note - in bash, when you want to put some of your code on the next line, you need to add a space and a 
# backslash ( \) to signal that the code continues onto the next line.

# <> <> <> <> <> #

# defaults only
WDIR=~/bigdata/250925_testrun
AmQUB_part1.sh -f "${WDIR}/data/dataAB.fastq" -p "${WDIR}/data/dataAB_map.txt"
AmQUB_part1.sh -f "${WDIR}/data/dataXY.fastq" -p "${WDIR}/data/dataXY_map.txt"
AmQUB_part2.sh -f "${WDIR}/part1_dataAB_output" -p "${WDIR}/data/dataAB_map.txt" 
AmQUB_part2.sh -f "${WDIR}/part1_dataAB_output" -p "${WDIR}/data/dataXY_map.txt" 

# defaults UPARSE
AmQUB_part3.sh --in "${WDIR}/dataAB_output"+"${WDIR}/dataXY_output" --len 300 \
	--out output_dir --al UPARSE

# defaults UNOISE3
AmQUB_part3.sh --in ID1_raw_output+ID2_raw_output --len 300 \
	--out output_dir --al UNOISE3

AmQUB_part4.sh -i output_dir -t 64 -e email@email.com

AmQUB_part5.sh -i output_dir

# <> <> <> <> <> #

# using options
WDIR=~/bigdata/250925_testrun
AmQUB_part1.sh -f "${WDIR}/data/AB.fastq" -p "${WDIR}/data/AB_map.txt" -m 2
AmQUB_part1.sh -f "${WDIR}/data/XY.fastq" -p "${WDIR}/data/XY_map.txt" -m 2
AmQUB_part2.sh -f "${WDIR}/part1_AB_output" -p "${WDIR}/data/AB_males_map.txt" -m 2 -r 200-301 -i 10 -s males
AmQUB_part2.sh -f "${WDIR}/part1_XY_output" -p "${WDIR}/data/XY_males_map.txt" -m 2 -r 200-301 -i 10 -s males

# all options UPARSE
AmQUB_part3.sh --in "part2_AB_males_output"+"part2_XY_males_output" --len 300 \
	--out output_dir_males --al UPARSE --fr yes --fe yes --fm 0.5 \
	--fd 2 --ug yes --us 2 --min 8 --ag yes --pid 0.98 \
	--map merged_map.txt --col Soil_Type+Tissue \
	--pre IDx_OTUs.fa --tblid 0.98 --un yes --rmf yes  

# all options UNOISE3
AmQUB_part3.sh --in "part2_AB_males_output"+"part2_XY_males_output" --len 300 \
	--out output_dir_males --al UPARSE --fr yes --fe yes --fm 0.5 \
	--fd 2 --ug yes --us 2 --min 8 --ag yes --alpha 2.5 \
	--map merged_map.txt --col Soil_Type+Tissue \
	--pre IDx_OTUs.fa --tblid 0.98 --un yes --rmf yes  

#    AmQUB_part4.sh  -i output_dir -t 256 -e email@email.com \
#        -b mega-blast -v 0.005 -c classifier.qza -f 0.8 --nostrategy1 --strategy2 --strategy3 
#    AmQUB_part4.sh  -o output_dir -t 256 -e email@email.com \
#        -c classifier.qza -f disable --strategy2 --strategy3 -s \ 

# <> <> <> <> <> #

# parameter files
WDIR=~/bigdata/250925_testrun
AmQUB_part1.sh "${WDIR}/params/part1_params.csv"
AmQUB_part2.sh "${WDIR}/params/part2_params.csv"
AmQUB_part3.sh "${WDIR}/params/part3_params.csv"
AmQUB_part4.sh "${WDIR}/params/part4_params.csv"
AmQUB_part5.sh "${WDIR}/params/part5_params.csv"
