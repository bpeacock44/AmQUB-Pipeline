#
# A small library of biom-format helper functions for QIIME work
# by Paul Ruegger
#
# To use:
# source qiime_helper_functions.sh
#

echo -e '* Type "biom_help" to see a list of biom-format helper functions *'

function biom_help {
    #print function names in this script
    echo
    echo -e '## BIOM functions that work in QIIME >= 1.9.1 ##'
    
    echo -e '# Convert biom (HDF5 or JSON) to txt:  biom2txt <in.biom> <out.txt>'
    echo -e '  Optionally (default=taxonomy):       biom2txt <in.biom> <out.txt> <header-key>'
    echo -e '# Convert (HDF5 or JSON) NO taxonomy:  biom2txt_notax <in.biom> <out.txt>'
    echo -e '# Convert text to biom (HDF5):         txt2biom <in.txt> <out.biom>'
    echo -e '  Optionally (default=taxonomy):       txt2biom <in.txt> <out.biom> <taxonomy|naive|sc_separated>'
    echo -e '# Convert text to biom (JSON):         txt2JSONbiom <in.txt> <out.biom>'
    echo -e '  Optionally (default=taxonomy):       txt2JSONbiom <in.txt> <out.biom> <taxonomy|naive|sc_separated>'
    
    echo -e '# Convert text to biom (HDF5) NO tax:  txt2biom_notax <in.biom> <out.txt>'
    echo -e '# Convert text to biom (JSON) NO tax:  txt2JSONbiom_notax <in.biom> <out.txt>'
    
    echo -e '# Convert biom (HDF5) to JSON:         biom2json <in.biom> <out.biom>'
    echo -e '# Convert biom (JSON) to HDF5:         json2biom <in.biom> <out.biom>'
    echo -e '# Add taxonomy from assignments:       biomAddObservations <in.biom> <out.biom> <tax_assignments_with_headers_file>'
    echo -e '# Summarize biom:                      biomsummarize <infile> <outfile>'
}


#Functions that work in QIIME >= 1.9.1:

function biom2txt {
    #$1 is the input biom file, $2 is the output file, $3 is a header in $1 (default="taxonomy") (optional)
    _K_="$3";
    
    
    #if $1 or $2 are NOT defined... notify user
    if [[ -z $1 ]] || [[ -z $2 ]]; then
        echo "Usage: biom2txt <table.biom> <table.txt> [<metadata>]";
        echo " <metadata> tag is optional (default=taxonomy)"
        kill -INT $$
    fi
    if ! [[ -e $1 ]]; then
        echo "Input error: file [$1] does not exist"'!';
        kill -INT $$
    fi
    
    #check filetype of biom
    FTYPE=$(file "$1" | perl -ne 'if(/Hierarchical Data Format/){print "hdf5"}elsif(/text.+long\slines/){print "json"}else{print ""}')
    #define a (portion of a) regular expression based on the filetype (helps later on to accurately determine if meta-data exists)
    if [[ $FTYPE == "hdf5" ]]; then
        RGX=""
    elif [[ $FTYPE == "json" ]]; then
        RGX=": {\""
    else
        echo "Input error: file [$1] doesn't look like an HDF5 or JSON biom file"'!'
        kill -INT $$
    fi
    
    #if $3 is NOT defined...
    if [[ -z $_K_ ]]; then
        #if here, $1 and $2 must be defined, so we will assume that taxonomy should be included in output
        #(but we'll check first if it's actually in $1 or not)
        if [[ $(grep -c "${RGX}taxonomy" $1) -gt 0 ]]; then
            #taxonomy was found
            echo "Meta-data header [taxonomy] found in biom. Saving with meta-data"
            biom convert -i "$1" -o "$2" --table-type="ASV table" --to-tsv --header-key "taxonomy";
        else
            #taxonomy was not found, so we won't try to include it (or we'll get a hanging 'taxonomy' header)
            echo "Meta-data header [taxonomy] was not found in biom. Saving without meta-data"
            biom convert -i "$1" -o "$2" --table-type="ASV table" --to-tsv;
        fi
    else #$1, $2 and $3 ARE defined
        #check if $3 is in $1
        if [[ $(grep -c "${RGX}$_K_" $1) -gt 0 ]]; then
            echo "Meta-data header [$_K_] found in biom. Saving with meta-data"
            biom convert -i "$1" -o "$2" --table-type="ASV table" --to-tsv --header-key "$_K_";
        else
            #notify but DO NOT save if user-defined $3 was NOT found (unlike case of meta-data default header, "taxonomy")
            echo -e "Meta-data header [$_K_] NOT FOUND in biom.\nFile [$2] NOT saved"'!';
            kill -INT $$
        fi
    fi
    
    #TODO: the logic below is a little screwed up. can't tell if $2 was JUST saved or already existed. Better to check at the outset?
    
    #if $1 AND $2 are defined, and $2 exists as a file, and '.0's are found in $2
    if ! [[ -z $1 ]] && ! [[ -z $2 ]] && [[ -e $2 ]] && [[ "$(grep -c "\.0\b" $2)" -gt 0 ]]; then
        #remove '.0' in the text file $2
        perl -ne 's/\.0\b//g; print;' < "$2" > "${2}a"; mv -f "${2}a" "${2}";
    fi
    #double-check $2 exists and notify
    if [[ -e $2 ]]; then echo "Saved [$2]"; else echo "File [$2] NOT saved"'!'" [header-key=$_K_]"; fi
    unset _K_ FTYPE RGX;
}

function biom2txt_notax {
    #if $1 or $2 are NOT defined... notify user
    if [[ -z $1 ]] || [[ -z $2 ]]; then
        echo "Usage: biom2txt_notax <table.biom> <table.txt>";
        echo " This function will not include any existing taxonomy"
        kill -INT $$
    fi
    if ! [[ -e $1 ]]; then
        echo "Input error: file [$1] does not exist"'!';
        kill -INT $$
    fi
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-tsv;
    perl -ne 's/\.0\b//g; print;' < "$2" > "${2}a"
    mv -f "${2}a" "${2}"
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}

function txt2biom {
    _K_="$3";#valid choices: [taxonomy|naive|sc_separated]
    if [ -z "$_K_" ]; then _K_='taxonomy'; else echo "setting header-key=$_K_"; fi
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-hdf5 --process-obs-metadata "$_K_";
    if [ -e "$2" ]; then echo "Saved [$2] [header-key=$_K_]"; else echo "File [$2] not saved"'!'" [header-key=$_K_]"; fi
    unset _K_;
}

function txt2biom_notax {
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-hdf5;
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}
function txt2JSONbiom_notax {
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-json;
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}

function txt2JSONbiom {
    _K_="$3";#valid choices: [taxonomy|naive|sc_separated]
    if [ -z "$_K_" ]; then _K_='taxonomy'; else echo "setting header-key=$_K_"; fi
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-json --process-obs-metadata "$_K_";
    if [ -e "$2" ]; then echo "Saved [$2] [header-key=$_K_]"; else echo "File [$2] not saved"'!'" [header-key=$_K_]"; fi
    unset _K_;
}

function biom2json {
    #if $1 or $2 are NOT defined... notify user
    if [[ -z $1 ]] || [[ -z $2 ]]; then
        echo "Usage: biom2json <table.biom> <table.json>";
        kill -INT $$
    fi
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-json;
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}

function json2biom {
    biom convert -i "$1" -o "$2" --table-type="ASV table" --to-hdf5;
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}

function biomsummarize {
    biom summarize-table -i "$1" -o "$2"
    if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
}

function biomAddObservations {
    #if $1, $2 or $3 are NOT defined...
    if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
        echo "Usage: biomAddObservations <in.biom> <out.biom> <meta-data-file>";
        echo "(Note that you must add headers to the meta data file if there are none. Eg., '#ASV ID<tab>taxonomy<tab>confidence')"
    else
        biom add-metadata -i "$1" -o "$2" --observation-metadata-fp "$3" --sc-separated taxonomy
        if [ -e "$2" ]; then echo "Saved [$2]"; else echo "File [$2] not saved"'!'; fi
    fi
}

#######################################
# Function for beta_diversity cleanup #
function condense_multiple_PCOA_plots {  #qiime_shell_helper_functions.sh
#
#Each folder containing a QIIME-created PCOA plot contains a subfolder, "emperor_required_resources", that contains required plotting functions
#In cases where you have multiple PCOA plots, each will contain its own copy of the emperor_required_resources subfolder
#This function renames/moves all PCOA index.html files into a single folder, then deletes all redundant emperor_required_resources subfolders,
# +resulting in a saving of ~1.4Mb of disk space per deleted copy
#
    #if $1 IS defined we will prepend the output .html files with it
    if ! [[ -z $1 ]]; then
        _GROUP_="$1"
        echo "_GROUP_=[$_GROUP_]"
    fi
    #define subfolder vars
    PLOTTING_RESOURCES_DIR="emperor_required_resources"
    PCOA_DIR="./pcoa_plots"
    
    #ensure that at least one <betadiv_metric>_emperor_pcoa_plot directory exists in the current directory
    EMPEROR_PCOA_PLOT_DIRS=($(find . -mindepth 1 -maxdepth 1 -type d -name "*_emperor_pcoa_plot" -printf '%P\n'))
    
    if [[ ${#EMPEROR_PCOA_PLOT_DIRS[@]} -ge 1 ]]; then
        #make subfolder into which all PCOA plots will be moved
        mkdir -v $PCOA_DIR
        mkdir $PCOA_DIR/dms; #for dm, pc and log files
        
        #ensure that $PCOA_DIR was created, then copy renamed versions of index.html files into it
        if [[ -d $PCOA_DIR ]]; then
            echo "Subdir [$PCOA_DIR] created";
            #put a copy of the plotting resources into the collection subdir
            echo -e "Copying a [$PLOTTING_RESOURCES_DIR] directory to [$PCOA_DIR]\n"
            cp -r ${EMPEROR_PCOA_PLOT_DIRS[0]}/$PLOTTING_RESOURCES_DIR/ $PCOA_DIR/
            #set timestamp of copy to the original's
            touch $PCOA_DIR/* -r ${EMPEROR_PCOA_PLOT_DIRS[0]}/$PLOTTING_RESOURCES_DIR
            
            #cycle through <betadiv_metric>_emperor_pcoa_plot dirs, copying/renaming the index.html files
            for EMP_DIR in ${EMPEROR_PCOA_PLOT_DIRS[@]}; do
                # EMP_DIR=${EMPEROR_PCOA_PLOT_DIRS[0]}; echo $EMP_DIR
                EMP_DIR2=${_GROUP_}${EMP_DIR};
                echo "Copying [$EMP_DIR/index.html] to [$PCOA_DIR/${EMP_DIR2}.html]"
                cp $EMP_DIR/index.html $PCOA_DIR/${EMP_DIR2}.html
                touch $PCOA_DIR/${EMP_DIR2}.html -r $EMP_DIR/index.html
            done
            
            #find and rename ancillary files
            DM_PC_LOGS=($(find . -mindepth 1 -maxdepth 1 -type f -regex '.*[a-z]+.*\([mc]\|[0-9]\)\.txt' -printf '%P\n'))
            if [[ ${#DM_PC_LOGS[@]} -gt 0 ]]; then
                rename -v 's//'$_GROUP_'/' ${DM_PC_LOGS[@]}
                #refind
                DM_PC_LOGS=($(find . -mindepth 1 -maxdepth 1 -type f -regex ".*$_GROUP_.*\.txt" -printf '%P\n'))
                #move ancillary files
                for _F_ in ${DM_PC_LOGS[@]}; do
                    mv -vf ${_F_} $PCOA_DIR/dms/
                done
            fi
            
            #verify the html files were copied properly then remove originals (TODO: probably don't need -s check?)
            for EMP_DIR in ${EMPEROR_PCOA_PLOT_DIRS[@]}; do
                EMP_DIR2=${_GROUP_}${EMP_DIR};
                if [[ -e $PCOA_DIR/${EMP_DIR2}.html ]] && [[ -s $PCOA_DIR/${EMP_DIR2}.html ]]; then
                    echo "(Removing [$EMP_DIR/])"
                    rm -rf $EMP_DIR
                fi
            done
        else
            echo "*** Error *** Unable to create subdir [$PCOA_DIR]"'!'" :(";
        fi
    else
        echo "*** No <betadiv_metric>_emperor_pcoa_plot directories exist in the current directory ***"
        echo "Please cd to a beta_diversity subfolder before calling this function"
    fi
    unset PLOTTING_RESOURCES_DIR PCOA_DIR EMPEROR_PCOA_PLOT_DIRS EMP_DIR _GROUP_ EMP_DIR2
}
# Function for beta_diversity cleanup #
#######################################


function count_taxa_levels {
    if [[ -z "$1" ]]; then
        echo "Counts occurences of k__, p__, c__, etc. in a taxonomic assignments file"
        echo "Usage: count_taxa_levels <taxonomic_assignments_file>"
        kill -INT $$
    fi
    _LEVELS_=("(k|d)__" p__ c__ o__ f__ g__ s__ u__ Unassigned)
    echo -e "Level\tCount"
    for _LEV_ in ${_LEVELS_[@]}; do
        _C_=$(grep -cP "$_LEV_" $1)
        echo -e "$_LEV_\t$_C_"
    done
    unset _LEVELS_ _LEV_ _C_
}



##################################################
# Print function names when this file is sourced #
function print_functions {
    # Find the full path to this script [from http://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself?rq=1]
    _PATH2QIIMEHELPERSCRIPT_=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)/`basename "${BASH_SOURCE[0]}"`
    # Get a list of the functions
    _MYFUNCIONS_=( $(grep -P "^function" ${_PATH2QIIMEHELPERSCRIPT_} | grep -vP "(grep|_stopUnless)" | awk '{print $2}') )
    # Format them with surrounding brackets
    _X_=$(printf "[%s]  " "${_MYFUNCIONS_[@]}");
    # Print them
    echo -e "Helper functions list:\n$_X_";
    # Cleanup the variables used here
    unset _MYFUNCIONS_ _X_ _PATH2QIIMEHELPERSCRIPT_ _MYFUNCIONS_
}
print_functions
# Print function names when this file is sourced #
##################################################