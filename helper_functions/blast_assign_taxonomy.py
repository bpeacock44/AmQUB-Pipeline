#!/usr/bin/env python3


"""
blast_assign_taxonomy.py

Loads an [ASVs2filter.log] file created by [blast_taxa_categorizer.py] and assigns taxonomy to the ASVs/query sequences

Usage:
    blast_assign_taxonomy.py -i <file> --db <file> -m <email@email.com> [-o <file>] [-K <int>] [-E <int>] [--add_sizes] [--assign_all] [--regex <regex>] [--asv <asv>]
    blast_assign_taxonomy.py --version
    blast_assign_taxonomy.py -h

Arguments:
    -i <file>     filepath to ASVs2filter.log
    --db <file>   filepath to json taxonomy database (eg: path-to/taxonomyDB.json)
    -m <email>    user's email for NCBI reporting

Options:
    -o <file>        output filepath (if not defined, prints to STDOUT)
    -K <int>         minimum acceptable bitscore for Keeper hits [default: 0]
    -E <int>         minimum acceptable bitscore for Environmental hits [default: 0]
    --assign_all     assign taxonomy even to Rejected ASVs (including PhiX) [default: False]
    --add_sizes      add ASV sizes to output [default: False]
    --regex <regex>  use a custom regular expression to detect an "Asv" line [default: ^Asv]
    --version        version
    --asv <asv>      debugging info
    -h               this help

Related scripts:
1) blast_taxa_categorizer.py (creates the [ASVs2filter.log] file)
"""


from collections import defaultdict
from urllib.error import HTTPError
from docopt import docopt
from Bio import Entrez
import xmltodict
import pickle
import json
import time
import sys
import re
import os


class ASV:
    '''Holds all blastn hits to a single asv found in an [ASVs2filter.log] file'''
    def __init__(self, asv, size, designation, predesignation):
        self.id = str(asv)
        self.size = int(size)
        self.designation = designation       # K or R
        self.predesignation = predesignation # EKU, N, etc.
        self.hits = {}
        self.best_hit_taxIDs = {}
        self.best_hit_bitscore = {}
        self.best_hit_pident = {}
        self.best_hit_qcov = {}
    def id(self):
        return self.id
    def size(self):
        return self.size
    def designation(self):
        return self.designation
    def predesignation(self):
        return self.predesignation
    def add_blast_hits(self, hitType, hitList):
        self.hits[hitType] = hitList
        self._set_best_hits(hitType, hitList)
    def _set_best_hits(self, hitType, hitList):
        best_taxIDs = []
        bestbitscore = 0.0
        best_pident = 0.0
        best_qcov = 0.0

        for hit in hitList:
            #eg. hit="111|87.879|1766|4|100"
            bitscore,pident,taxID,nHits,qcov = re.split("\|", hit)
            bitscore = float(bitscore)
            if bitscore > bestbitscore:
                bestbitscore = bitscore
                best_pident = float(pident)
                best_qcov = float(qcov)
        for hit in hitList:
            bitscore,pident,taxID,nHits,qcov = re.split("\|", hit)
            bitscore = float(bitscore)
            if bitscore == bestbitscore:
                best_taxIDs.append(taxID)
            else:
                break
        self.best_hit_taxIDs[hitType] = best_taxIDs
        self.best_hit_bitscore[hitType] = bestbitscore
        self.best_hit_pident[hitType] = best_pident
        self.best_hit_qcov[hitType] = best_qcov
    def get_hits(self, hitType):
        if hitType in self.hits:
            return self.hits[hitType]
    def get_best_hit_taxIDs(self, hitType):
        if hitType in self.best_hit_taxIDs:
            return self.best_hit_taxIDs[hitType]
    def get_best_hit_bitscore(self, hitType):
        if hitType in self.best_hit_bitscore:
            return self.best_hit_bitscore[hitType]
    def get_best_hit_pident(self, hitType):
        if hitType in self.best_hit_pident:
            return self.best_hit_pident[hitType]
    def get_best_hit_qcov(self, hitType):
        if hitType in self.best_hit_qcov:
            return self.best_hit_qcov[hitType]
#

def get_and_check_opts(args):
    opts = docopt(args, version="blast_assign_taxonomy.py, version 0.4")

    #verify that input files exists
    optionsFail = False
    if not os.path.isfile(opts['-i']):
        print ("** Error ** Input file ["+opts['-i']+"] does not exist!", file=sys.stderr)
        optionsFail = True

    if os.path.isfile(opts['--db']):
        print("** Using db ["+opts['--db']+"] **", file=sys.stderr)
        pass
    else:
        print("taxonomy db will be stored as: ["+opts['--db']+"]")
        #print("** Could not find db ["+opts['--db']+"] **", file=sys.stderr)
        #reply = str(input("Would you like to [c]reate a new db here or [q]uit?\n")).lower().strip()
        #if reply in ['c', 'q']:
        #    if reply in ['c']:
        #        pass
        #    else:
        #        print ("Quitting", file=sys.stderr)
        #        optionsFail = True
        #else:
        #    print("** Invalid choice. Please try again **", file=sys.stderr)
        #    optionsFail = True

    # Check for the email address argument
    if '-m' in opts and opts['-m']:
        email = opts['-m']
        # Validate email format using regular expression
        if not re.match(r"[^@]+@[^@]+\.[^@]+", email):
            print("** Error ** Invalid email address format. Please provide a valid email address using the -m option.", file=sys.stderr)
            optionsFail = True

    else:
        print("** Error ** Email address is required for posting to NCBI. Please provide it using the -m option.", file=sys.stderr)
        optionsFail = True

    if optionsFail:
        sys.exit(1)  # Exit with a non-zero status to indicate an error

    return opts
#

def get_starting_line(opts):
    starting_line = 0
    with open(opts['-i'], 'r') as infile:
        #skip irrelevant header info
        for line in infile:
            if re.search("->", line) is None:
                starting_line += 1
                continue
            else:
                break
    return (starting_line)


def parse_ASVs2filter_log(opts):
    ''' '''
    asv_list = []
    # stopAt = 5000 ###

    starting_line = get_starting_line(opts)

    with open(opts['-i'], 'r') as infile:

        #skip any irrelevant header info
        for _ in range(starting_line):
            next(infile)

        #-- real data processing starts here --#
        nASVs = 0
        finishedASV = False
        for line in infile:
            line = line.rstrip()
            if re.search(opts['--regex'], line) is not None:
                # Asv96762_2 [EKU]->[K]
                #create a new asvObj
                asvline = re.split(" ", line)                 # ['Asv96762_2', '[EKU]->[K]']
                asv,size = re.split("_", asvline[0])          # ['Asv96762_2', '2']
                designations = (re.split("->", asvline[1]))   # ['[EKU]', '[K]']
                predesignation = designations[0][1:-1]        # 'EKU'
                designation = designations[1][1:-1]           # 'K'
                asvObj = ASV(asv,size,designation,predesignation)
                continue
            elif re.search(r"^\S+size=", line) is not None:
                # F001_1027355;size=8;_0 [EKU]->[K]
                #create a new asvObj
                asvline = re.split(" ", line)                 # ['F001_100692;size=100;_0', '[EK]->[K]']
                asv,size = re.split(";_", asvline[0])         # ['F001_1027355;size=8', '0']
                asv = asv+";"
                designations = (re.split("->", asvline[1]))   # ['[EK]', '[K]']
                predesignation = designations[0][1:-1]        # 'EK'
                designation = designations[1][1:-1]           # 'K'
                asvObj = ASV(asv,size,designation,predesignation)
                continue
            if re.search("^TRINITY", line) is not None:
                #TRINITY_DN393_c0_g1_i1_27077 [KU]->[K]
                #create a new asvObj
                asvline = re.sub(r'_(\d+)', r'__\1', line)    #  'TRINITY_DN393_c0_g1_i1__27077 [KU]->[K]'
                asvline = re.split(" ", asvline)              # ['TRINITY_DN393_c0_g1_i1__27077', '[KU]->[K]']
                asv,size = re.split("__", asvline[0])         # ['TRINITY_DN393_c0_g1_i1', '27077']
                designations = (re.split("->", asvline[1]))   # ['[KU]', '[K]']
                predesignation = designations[0][1:-1]        # 'KU'
                designation = designations[1][1:-1]           # 'K'
                asvObj = ASV(asv,size,designation,predesignation)
                continue
            elif re.search(r"  E\[", line) is not None:
                hitType = "E"
            elif re.search(r"  K\[", line) is not None:
                hitType = "K"
            elif re.search(r"  R\[", line) is not None:
                hitType = "R"
            elif re.search(r"  U\[", line) is not None:
                hitType = "U"
                finishedASV = True
            else:
                continue

            #get line ready to add to asvObj
            line = line[4:]
            line = re.sub(r'[()\[\]]', '', line)
            hitList = re.split(" ", line)

            #add hitType and hitList to asvObj
            if len(line) > 0 and hitType is not None:
                asvObj.add_blast_hits(hitType, hitList)
                hitType = None

            if finishedASV:
                asv_list.append(asvObj)
                nASVs = nASVs + 1
                finishedASV = False
                # if nASVs >= stopAt:
                    # break

        return asv_list
#

def parse_xml_file(taxa_xml_file):
    #@Returns: dict[TaxId]=taxonomy
    d = {}

    # debug_mode = True
    debug_mode = False
    case1 = 0 ##debug
    case2 = 0 ##debug
    case3 = 0 ##debug
    case4 = 0 ##debug
    case5 = 0 ##debug

    #create a lookup dict for ranks
    ranks2keep = {}
    taxa_level_prefix = "k__ p__ c__ o__ f__ g__ s__ u__".split() #u__ = subspecies
    ranks = "superkingdom phylum class order family genus species subspecies".split()
    for rank in ranks:
        ranks2keep[rank] = 1

    #parse downloaded xml file to a dict
    with open(taxa_xml_file) as fh:
        rootD = xmltodict.parse(fh.read())

    if isinstance(rootD, dict) and isinstance(rootD['TaxaSet'], dict):
        taxon_objectsL = rootD['TaxaSet']['Taxon']
    else:
        print(":  ** Nothing to parse in [",taxa_xml_file,"] **", file=sys.stderr)
        return d

    #if there is only one taxon object in the XML file it gets returned as a dict instead of a list of dicts
    if isinstance(taxon_objectsL, dict):
        #so make it a list to keep downstream code happy
        taxon_objectsL = [taxon_objectsL]

    #
    for TaxonD in taxon_objectsL:
        TaxId = TaxonD['TaxId']                   #the eposted/efetched taxID
        ScientificName = TaxonD['ScientificName'] #the species/subspecies name
        taxonL = TaxonD['LineageEx']['Taxon']     #an ordered list of dicts, each containing rank info

        if isinstance(taxonL, dict):
            print(taxa_xml_file,": taxonL = [taxonL]")
            taxonL = [taxonL]

        taxonomy = []
        taxonomyD = {}
        genus = ""
        species = ""
        subspecies = ""
        species_found = False


        #initialize keys for desired taxonomic ranks
        for rank in ranks:
            if rank not in taxonomyD:
                taxonomyD[rank] = None

        #replace 'None' values if they are defined in the xml
        for taxlevD in taxonL:
            rank = taxlevD['Rank']
            if rank in taxonomyD:
                taxonomyD[rank] = taxlevD['ScientificName']

        #construct the full taxonomic lineage for the current Taxon in a greengenes-like format
        last_defined_taxa = ""
        for i,rank in enumerate(ranks):

            #use the genus name to help format a species/subspecies name from the ScientificName
            # (note that we'll only use this species name if a 'species' rank is not found)
            if rank == 'genus':
                if taxonomyD[rank] != None:
                    ## Example:
                    ## Genus = "Leisingera"
                    ## ScientificName = "Leisingera methylohalidivorans DSM 14336"

                    genus = taxonomyD[rank]

                    #remove the genus from the (species_subspecies) ScientificName
                    species_subspecies = re.sub(r""+genus+" ", "", ScientificName)
                    ## species_subspecies = "methylohalidivorans DSM 14336"

                    #split the remainder on whitespace
                    species_subspecies = re.split(" ", species_subspecies)

                    #the species name is the genus + the first word(?) in species_subspecies
                    species = genus + "_" + species_subspecies[0]

                    #everything else is the subspecies (AFAIK)
                    subspecies = "".join(species_subspecies[1:])
                    ## subspecies = "DSM14336"

                    #add to taxonomy dict (provisionally)
                    taxonomyD['species2'] = species
                    taxonomyD['subspecies2'] = subspecies
                else:
                    taxonomyD['species2'] = None
                    taxonomyD['subspecies2'] = None


            #if a species rank is found (sometimes it isn't) find a/the subspecies name
            if rank == 'species':
                if taxonomyD[rank] != None:
                    species_found = True
                    #this info will supercede any "species" found above (if rank == 'genus')
                    species = taxonomyD[rank]
                    #try to find a subspecies name by removing the species name from the ScientificName
                    subspecies = re.sub(r""+species+" ", "", ScientificName)
                    #remove any remaining whitespace
                    subspecies = re.sub(r" ", "", subspecies)
                    #add to taxonomy dict
                    taxonomyD['subspecies'] = subspecies
                else:
                    #use what we found above, if anything
                    taxonomyD['species']    = taxonomyD['species2']
                    taxonomyD['subspecies'] = taxonomyD['subspecies2']

            if rank == 'subspecies':
                if taxonomyD[rank] == None:
                    taxonomyD['subspecies'] = re.sub(r" ", "_", ScientificName)


            #append the current rank's value to the taxonomy list
            if taxonomyD[rank] is not None:
                taxonomyD_rank = re.sub(r" ", "_", taxonomyD[rank])
                taxonomy.append(taxa_level_prefix[i]+taxonomyD_rank)
                last_defined_taxa = taxonomyD_rank
            else:
                taxonomy.append(taxa_level_prefix[i]+"unclassified_"+last_defined_taxa)


        ### ### ###
        if debug_mode:
            #case1: species_found = True, species = methylohalidivorans, subspecies = DSM14336  #<--definitely happens
            #case2: species_found = True, taxonomyD['species'] != methylohalidivorans           #<--possibly happens
            #case3: species_found = False, species = methylohalidivorans, subspecies = DSM14336 #<--definitely happens
            #case4: species_found = False, species = methylohalidivorans, subspecies = ""       #<--definitely happens
            #case5: species_found = False, species = "", subspecies = ""                        #<--definitely happens

            #case1
            #after checking all existing Ranks
            if species_found and taxonomyD['species'] == species:
                case1 += 1

            #case2
            if species_found and taxonomyD['species'] != species:
                case2 += 1
                sys.stderr.write("case2: %s # %s\n" % (taxonomyD['species'],species))
                sys.stderr.write("%s\n" % taxonomy)
                sys.exit(0)

            #case3
            if not species_found and len(species) > 0 and len(subspecies) > 0:
                case3 += 1

            #case4
            if not species_found and len(species) > 0 and len(subspecies) == 0:
                case4 += 1

            #case5
            if not species_found and len(species) == 0 and len(subspecies) == 0:
                case5 += 1
                # if case5 == 2:
                sys.stderr.write("case5: %s # %s # %s #" % (TaxId,ScientificName,species))
                sys.stderr.write("%s\n" % taxonomy)
                    # sys.exit(0)
        ### ### ###


        #remove whitespace and add to dict
        lineage = ";".join(taxonomy).replace(" ","")
        lineage = lineage.replace("=",".")
        lineage = lineage.replace("'","")

        #add to dict
        d[TaxId] = lineage

        #if any AkaTaxIds exist in taxa_xml_file, add them too
        if 'AkaTaxIds' in TaxonD:
            AkaTaxIds = re.split(",", TaxonD['AkaTaxIds']['TaxId'])
            for AkaTaxId in AkaTaxIds:
                d[AkaTaxId] = lineage
        """
        #TODO: check whether multiple AkaTaxIds are separated by a "," or something else
        #From the taxa_xml_file:
        "AkaTaxIds": {
            "TaxId": "1698445" #<-- ie., how will these be separated if there are more than one?
        },
        """
    #
    if debug_mode:
        sys.stderr.write("case1: %s  case:2 %s  case3: %s  case4: %s  case5: %s\n" % (case1,case2,case3,case4,case5))

    return d
#

def get_localDB_taxonomy(opts):
    #@Returns dict[taxonIDs]=taxonomic_lineages
    old_taxid_dict = {}
    if os.path.isfile(opts['--db']):
        with open(opts['--db'], 'r') as fh:
            old_taxid_dict = json.load(fh)

    return old_taxid_dict
#

def get_ASVs2filter_taxonIDs_list(opts, asv_list):
    #@Returns list[taxonIDs]
    taxIDs = []
    for asv in asv_list:
        #we'll always want "Keeper" taxIDs
        if asv.designation == "K":
            if asv.get_best_hit_taxIDs("K") is not None:
                taxIDs += asv.get_best_hit_taxIDs("K")
            if asv.get_best_hit_taxIDs("E") is not None:
                taxIDs += asv.get_best_hit_taxIDs("E")

        #and sometimes we'll want all taxIDs
        if opts["--assign_all"]:
            #we may assign from R, RU, EKRU, ... predesignations
            if asv.get_best_hit_taxIDs("R") is not None:
                taxIDs += asv.get_best_hit_taxIDs("R")
            #ensure we include "Undesignated" taxa
            if "U" in asv.predesignation:
                taxIDs += asv.get_best_hit_taxIDs("U")

    taxIDs_dict = {}
    for taxID in taxIDs:
        taxIDs_dict[taxID] = 1

    return list(taxIDs_dict.keys())
#

def epost_taxonIDs(opts, new_taxonIDs_list):
    #TODO: define these earlier? somewhere else?
    #set some options
    opts['taxa_xml_file'] = "taxa_01." # TODO: define directory as /tmp
    opts['retmax'] = 1000
    opts['retmode'] = "xml"
    opts['max_attempts'] = 3
    opts['count'] = len(new_taxonIDs_list)

    if opts['count'] > 0:
        print(":  ",str(opts['count'])+" taxIDs will be ePost'ed to NCBI", file=sys.stderr)
        #epost taxon ID list
        opts['db'] = "taxonomy";
        Entrez.email = opts['email']
        epost_results = Entrez.read(Entrez.epost(opts['db'], id=",".join(new_taxonIDs_list)))
        #get WebEnv and QueryKey from epost_results
        opts['WebEnv']   = epost_results["WebEnv"]
        opts['QueryKey'] = epost_results["QueryKey"]
#

def get_next_xml_idx(opts):
    #finds and returns an index 1 greater than any existing xml files
    idx = 1
    xml_file_exists = os.path.isfile(opts['taxa_xml_file']+str(idx)+".xml")
    #find the last one that exists
    while xml_file_exists:
        idx += 1
        xml_file_exists = os.path.isfile(opts['taxa_xml_file']+str(idx)+".xml")

    return idx
#

def download_eposted_taxonIDs_to_XML(opts):
    # download taxonomies in XML format
    opts['xml_files'] = []
    idx = get_next_xml_idx(opts)
    for start in range(0, opts['count'], opts['retmax']):
        out_file = opts['taxa_xml_file'] + str(idx) + ".xml"
        with open(out_file, 'w') as out_handle:
            end = min(opts['count'], start + opts['retmax'])
            num2fetch = str(end - start)
            attempt = 1
            while attempt <= opts['max_attempts']:
                print(": Efetching taxa.")
                try:
                    #print(":  Efetching [" + num2fetch + "] taxa. Attempt", attempt, file=sys.stderr)
                    fetch_handle = Entrez.efetch(db=opts['db'], rettype="", retmode=opts['retmode'],
                                                 retstart=start, retmax=opts['retmax'],
                                                 webenv=opts['WebEnv'], query_key=opts['QueryKey'])
                    data = fetch_handle.read()
                    fetch_handle.close()
                    out_handle.write(data.decode('utf-8'))
                    opts['xml_files'].append(out_file)
                    idx += 1
                    break
                except HTTPError as err:
                    if 500 <= err.code <= 599 or err.code == 400:
                        print(":  Received error from server %s" % err, file=sys.stderr)
                        attempt += 1
                        time.sleep(15 * attempt)  # Increasing delay with each attempt
                        continue
                    else:
                        print(":  Error from server %s" % err, file=sys.stderr)
                        raise

#

def get_taxonomies_from_XML_files(opts):
    #@Returns dict[taxonIDs]=taxonomic_lineages
    all_new_taxonomies_dict = {}
    for out_file in opts['xml_files']:
        #parse the downloaded XML file into a dict
        print(":  Parsing XML file ["+out_file+"]", file=sys.stderr)
        new_taxIDs_dict = parse_xml_file(out_file)
        #print(":  Parsed",len(new_taxIDs_dict),"new TaxIDs", file=sys.stderr)

        #merge the old and new dicts
        #print(":    (pre-update) len(all_new_taxonomies_dict:",len(all_new_taxonomies_dict), file=sys.stderr)
        all_new_taxonomies_dict.update(new_taxIDs_dict)
        #print(":    (post-update)len(all_new_taxonomies_dict:",len(all_new_taxonomies_dict), file=sys.stderr)

    return all_new_taxonomies_dict
#

def taxID_is_integer(x):
    non_num = re.search(r"\D+", x)
    if non_num:
        return False
    else:
        return True

#a control function to get taxonomies from all sources
def get_taxonomy(opts, asv_list):
    #@Returns dict[taxonIDs]=taxonomic_lineages

    localDB_taxonomy_dict = get_localDB_taxonomy(opts)
    ASVs2filter_taxonIDs_list = get_ASVs2filter_taxonIDs_list(opts, asv_list)

    #determine which of the ASVs2filter_taxonIDs are new
    new_taxonIDs_list = list(set(ASVs2filter_taxonIDs_list) - set(localDB_taxonomy_dict.keys()))

    #ensure all taxIDs are okay (are integers)
    new_taxonIDs_list[:] = [x for x in new_taxonIDs_list if taxID_is_integer(x)]

    if opts["debug"]:
        print(": len(localDB_taxonomy_dict):",len(localDB_taxonomy_dict), file=sys.stderr)
        print(": len(ASVs2filter_taxonIDs_list):",len(ASVs2filter_taxonIDs_list), file=sys.stderr)
        print(": len(new_taxonIDs_list):",len(new_taxonIDs_list), file=sys.stderr)

    #update localDB_taxonomy_dict with any new taxa
    if len(new_taxonIDs_list) > 0:
        epost_taxonIDs(opts, new_taxonIDs_list)
        download_eposted_taxonIDs_to_XML(opts)
        all_new_taxonomies_dict = get_taxonomies_from_XML_files(opts)
        localDB_taxonomy_dict.update(all_new_taxonomies_dict)
        with open(opts['--db'], 'w') as fh:
            json.dump(localDB_taxonomy_dict, fh, sort_keys=True, indent=4)

        if opts["debug"]:
            print(": len(localDB_taxonomy_dict):",len(localDB_taxonomy_dict), file=sys.stderr)
    else:
        if opts["debug"]:
            print(": No new taxIDs encountered in input file... using saved taxonomic info only", file=sys.stderr)

    #remove subspecies?
    if opts["--del_ssp"]:
        if opts["debug"]:
            print(": opts[--del_ssp]", file=sys.stderr)
        for k in list(localDB_taxonomy_dict.keys()):
            v = re.sub(r";u_\S+","", localDB_taxonomy_dict[k])
            localDB_taxonomy_dict[k] = v

    return localDB_taxonomy_dict
#

def devel__print_blast_hit_summaries(opts, asv_list, taxonomy_dict):
    ignored_taxonomy_dict = {}
    ignored_taxonomy_dict['77133'] = 1
    ignored_taxonomy_dict['378812'] = 1
    ignored_taxonomy_dict['258848'] = 1

    Kbits = {}
    Kbits2 = {}
    #add info to Kbits dict
    for asv in asv_list:
        #print(asv.id+"_"+str(asv.size), asv.designation)
        if asv.designation == "K":
            asv_name = asv.id+"_"+str(asv.size)
            if asv.get_best_hit_taxIDs("K") is not None:
                best_hit_bitscoreK = asv.get_best_hit_bitscore("K")
                #print("  Kbit["+str(best_hit_bitscoreK)+"]")
                #make a list of top taxonomies
                tmptaxa = []
                tmptaxa2 = []
                for txid in asv.get_best_hit_taxIDs("K"):
                    if txid not in ignored_taxonomy_dict:
                        #print("   "+txid+"\t"+taxonomy_dict[txid])
                        tmptaxa.append([str(txid)+"\t"+taxonomy_dict[txid]])
                        tmptaxa2.append(taxonomy_dict[txid].split(";"))
                        # tmptaxa2.append(taxonomy_dict[txid])

                #if there are >= N taxonomies for this asv
                if len(tmptaxa2) >= 2:
                    #add it tmptaxa to Kbits
                    if best_hit_bitscoreK not in Kbits:
                        Kbits[best_hit_bitscoreK] = {}
                        Kbits2[best_hit_bitscoreK] = {}
                    Kbits[best_hit_bitscoreK][asv_name] = tmptaxa
                    Kbits2[best_hit_bitscoreK][asv_name] = tmptaxa2
    #
    if 0:
        #prints taxonomy in Kbits
        print(len(list(Kbits.keys())),"bitscores found")
        bitscoresK = list(Kbits.keys())
        bitscoresK.sort(key=lambda x: float(x), reverse=True)
        for best_hit_bitscoreK in bitscoresK:
            print("Kbit["+str(best_hit_bitscoreK)+"]")
            for asv_name in Kbits[best_hit_bitscoreK]:
                print("  ["+str(asv_name)+"]")
                for taxa in Kbits[best_hit_bitscoreK][asv_name]:
                    print("    ",''.join(taxa))
            print()
    else:
        #finds lowest common taxonomic level
        bitscoresK = list(Kbits2.keys())
        bitscoresK.sort(key=lambda x: float(x), reverse=True)
        for best_hit_bitscoreK in bitscoresK:
            lowComTaxaLev = {}
            print("Kbit["+str(best_hit_bitscoreK)+"]")
            for asv_name in Kbits2[best_hit_bitscoreK]:
                # print("  ["+str(asv_name)+"]")
                # for taxa in Kbits2[best_hit_bitscoreK][asv_name]:
                    # print("    ",' '.join(taxa))
                taxonomies_list = Kbits2[best_hit_bitscoreK][asv_name]
                lowest_taxa = get_lowest_common_taxonomy2(opts, taxonomies_list, asv)
                lentxa = len(lowest_taxa)
                if lentxa not in lowComTaxaLev:
                    lowComTaxaLev[lentxa] = 1
                else:
                    lowComTaxaLev[lentxa] += 1
            lengths = list(lowComTaxaLev.keys())
            lengths.sort(key=lambda x: x, reverse=False)
            print("len(Kbits2[best_hit_bitscoreK]) =",len(Kbits2[best_hit_bitscoreK]))
            for l in lengths:
                print(l, lowComTaxaLev[l])

            print()
#

# def get_lowest_common_taxonomy(opts, taxonomies_list, asv):
# #Note: does not currently exclude "unassigned_..." taxonomies, which can be problematic
# #TODO try to remove assignments with 'unassigned_' in them (unless that's all we have)?
    # if len(taxonomies_list) == 0:
        # return []

    # # ####################
    # # ### experimental ###
    # # # remove unclassified_ assignments before finding LCT
    # # to_remove = []
    # # for i in range(0, len(taxonomies_list)):
        # # if re.search("g__unclassified", ';'.join(taxonomies_list[i])) is not None:
            # # to_remove.append(i)
    # # if len(to_remove) > 0 and len(to_remove) < len(taxonomies_list):
        # # for j in reversed(to_remove):
            # # taxonomies_list.pop(j)
    # # ### experimental ###
    # # ####################

    # #start with the first one in the list
    # lowest_common_taxa = taxonomies_list[0]

    # #return immediately if we only have 1 taxonomy
    # if len(taxonomies_list) == 1:
        # return lowest_common_taxa

    # #find lowest common taxonomy in list
    # for k in range(1, len(taxonomies_list)):
        # #lowest_common_taxa = [i for i, j in zip(lowest_common_taxa, taxonomies_list[k]) if i == j]
        # #but ignore a mismatch if it's due to an "__unclassified" taxonomic level (before the genus level)
        # # lowest_common_taxa = [i for i, j in zip(lowest_common_taxa, taxonomies_list[k]) if i == j or re.search("__unclassified", j) is not None] # !!: k__Bacteria;g__unclassified_Gemmatimonadaceae
        # # lowest_common_taxa = [ i for i, j in zip(lowest_common_taxa, taxonomies_list[k])
            # # if i == j and not (re.search("__unclassified", i) is not None and re.search("__unclassified", j) is not None)]
        # lowest_common_taxa = [ i for i, j in zip(lowest_common_taxa, taxonomies_list[k])
            # if i == j and (re.search("__unclassified", i) is None and re.search("__unclassified", i) is None) ]

    # ###pdb
    # if asv.id == "Asv701":
        # lct = lowest_common_taxa
        # print("asv.id =", asv.id)
    # ###pdb
    # ##no common taxonomic level was found! try to find a majority consensus##
    # if len(lowest_common_taxa) == 0:
        # dc = defaultdict(lambda: defaultdict(lambda: 0))
        # #get min length of all lists inside taxonomies_list
        # minrange = min(map(len, taxonomies_list)) + 1
        # levels = [0]
        # for level in range(1, minrange):
            # levels.append(level)
            # # if levels[-1] == 0:
                # # continue
            # #print("#", level, levels, levels[-1], "#")
            # for txa in taxonomies_list:
                # cumTaxa = ';'.join(txa[0:levels[-1]]);
                # dc[level][cumTaxa] += 1
                # #print("##", level, "|", cumTaxa, "|", dc[level][cumTaxa], "##")

        # #decide
        # for level in range(1, minrange):
            # #get taxonomy
            # maxTaxakey = max(dc[level], key=dc[level].get)
            # maxTaxaCount = dc[level][maxTaxakey]
            # maxTaxaCountCount = sum(value == maxTaxaCount for value in dc[level].values())
            # #print(level, maxTaxakey, maxTaxaCount, maxTaxaCountCount)
            # if maxTaxaCountCount > 1:
                # if level == 1:
                    # #disagreement at k__
                    # #print("Return: k__Undetermined", maxTaxaCountCount, "#1")
                    # #lowest_common_taxa = ["k__Undetermined"]
                    # allTaxaAtMaxCount = [taxa for taxa in dc[level].keys() if dc[level][taxa] == maxTaxaCount]
                    # allTaxa = "_OR_".join(allTaxaAtMaxCount)
                    # #print("Return:", allTaxa, "#1")
                    # lowest_common_taxa = [allTaxa]
                # else:
                    # #disagreement at p__ or lower, so we'll report the previous level
                    # backOneLevel = level - 1
                    # maxTaxakey = max(dc[backOneLevel], key=dc[backOneLevel].get)
                    # maxTaxaCount = dc[backOneLevel][maxTaxakey]
                    # maxTaxaCountCount = sum(value == maxTaxaCount for value in dc[backOneLevel].values())
                    # #maxTaxaCountTaxa = [taxa for taxa in dc[backOneLevel].keys() if dc[backOneLevel][taxa] == maxTaxaCount]
                    # #print("Return:", backOneLevel, "|", maxTaxakey, "|", maxTaxaCount, maxTaxaCountCount, "#2")
                    # if backOneLevel == 1:
                        # lowest_common_taxa = [maxTaxakey]
                    # else:
                        # lowest_common_taxa = ";".split(maxTaxakey)
                    # # print("halting at", asv.id)
                    # # pdb.set_trace()
                # break

    # # if(asv.id == opts["--asv"]):
        # # print("halting")
        # # pdb.set_trace()

    # return lowest_common_taxa
# #

def get_lowest_common_taxonomy2(opts, taxonomies_list, asv):
#Note: does not currently exclude "unassigned_..." taxonomies, which can be problematic
#TODO try to remove assignments with 'unassigned_' in them (unless that's all we have)?
    if len(taxonomies_list) == 0:
        return []

    # ####################
    # ### experimental ###
    # # remove unclassified_ assignments before finding LCT
    # to_remove = []
    # for i in range(0, len(taxonomies_list)):
        # if re.search("g__unclassified", ';'.join(taxonomies_list[i])) is not None:
            # to_remove.append(i)
    # if len(to_remove) > 0 and len(to_remove) < len(taxonomies_list):
        # for j in reversed(to_remove):
            # taxonomies_list.pop(j)
    # ### experimental ###
    # ####################

    #start with the first one in the list
    lowest_common_taxa = taxonomies_list[0]

    #return immediately if we only have 1 taxonomy
    if len(taxonomies_list) == 1:
        return lowest_common_taxa
    # #add a duplicate so we can detect and remove any "__unclassified" assignments??
    # if len(taxonomies_list) == 1:
        # taxonomies_list.append(taxonomies_list[0])

    #find lowest common taxonomy in list
    for k in range(1, len(taxonomies_list)):
        #lowest_common_taxa = [i for i, j in zip(lowest_common_taxa, taxonomies_list[k]) if i == j]
        #but ignore a mismatch if it's due to an "__unclassified" taxonomic level (before the genus level)
        # lowest_common_taxa = [i for i, j in zip(lowest_common_taxa, taxonomies_list[k]) if i == j or re.search("__unclassified", j) is not None] # !!: k__Bacteria;g__unclassified_Gemmatimonadaceae
        # lowest_common_taxa = [ i for i, j in zip(lowest_common_taxa, taxonomies_list[k])
            # if i == j and not (re.search("__unclassified", i) is not None and re.search("__unclassified", j) is not None) ]
        # lowest_common_taxa = [ i for i, j in zip(lowest_common_taxa, taxonomies_list[k])
            # if i == j and (re.search("__unclassified", i) is None and re.search("__unclassified", j) is None) ]
        lowest_common_taxa = [ i for i, j in zip(lowest_common_taxa, taxonomies_list[k])
            if i == j and (re.search("__unclassified", i) is None) ]

    ##no common taxonomic level was found! try to find a majority consensus##
    if len(lowest_common_taxa) == 0:
        dc = defaultdict(lambda: defaultdict(lambda: 0))
        #get min length of all lists inside taxonomies_list
        minrange = min(map(len, taxonomies_list)) + 1
        levels = [0]
        for level in range(1, minrange):
            levels.append(level)
            # if levels[-1] == 0:
                # continue
            #print("#", level, levels, levels[-1], "#")
            for txa in taxonomies_list:
                cumTaxa = ';'.join(txa[0:levels[-1]]);
                dc[level][cumTaxa] += 1

        #decide
        for level in range(1, len(taxonomies_list)):
            #get taxonomy
            maxTaxakey = max(dc[level], key=dc[level].get)
            maxTaxaCount = dc[level][maxTaxakey]
            maxTaxaCountCount = sum(value == maxTaxaCount for value in dc[level].values())
            #print(level, maxTaxakey, maxTaxaCount, maxTaxaCountCount)
            if maxTaxaCountCount >= 1:
                if level == 1:
                    #disagreement at k__
                    #print("Return: k__Undetermined", maxTaxaCountCount, "#1")
                    #lowest_common_taxa = ["k__Undetermined"]
                    #allTaxaAtMaxCount = [taxa for taxa in dc[level].keys() if dc[level][taxa] == maxTaxaCount]
                    #allTaxa = "_OR_".join(allTaxaAtMaxCount)
                    allTaxaAtMaxCount = [taxa for taxa in dc[level].keys()]
                    allTaxa = "_OR_".join(sorted(allTaxaAtMaxCount))
                    # print("Return:", allTaxa, maxTaxaCount, "#1")
                    # pdb.set_trace() ###pdb
                    lowest_common_taxa = [allTaxa]
                else:
                    #disagreement at p__ or lower, so we'll report the previous level
                    backOneLevel = level - 1
                    maxTaxakey = max(dc[backOneLevel], key=dc[backOneLevel].get)
                    maxTaxaCount = dc[backOneLevel][maxTaxakey]
                    maxTaxaCountCount = sum(value == maxTaxaCount for value in dc[backOneLevel].values())
                    #maxTaxaCountTaxa = [taxa for taxa in dc[backOneLevel].keys() if dc[backOneLevel][taxa] == maxTaxaCount]
                    print("Return:", backOneLevel, "|", maxTaxakey, "|", maxTaxaCount, maxTaxaCountCount, "#2")
                    if backOneLevel == 1:
                        lowest_common_taxa = [maxTaxakey]
                    else:
                        lowest_common_taxa = maxTaxakey.split(";")
                break

    return lowest_common_taxa
#

def get_consensus_taxonomy(opts, taxonomies_list):
    if len(taxonomies_list) == 0:
        return []

    #return immediately if we only have 1 taxonomy
    if len(taxonomies_list) == 1:
        return taxonomies_list[0]

    ####################
    ### experimental ###
    # remove unclassified_ assignments before finding consensus
    to_remove = []
    for i in range(0, len(taxonomies_list)):
        if re.search("g__unclassified", ';'.join(taxonomies_list[i])) is not None:
            to_remove.append(i)
    if len(to_remove) > 0 and len(to_remove) < len(taxonomies_list):
        for j in reversed(to_remove):
            taxonomies_list.pop(j)
    ### experimental ###
    ####################

    #"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Xanthomonas;s__Xanthomonas_euvesicatoria;u__pv.allii"
    #get counts of each genera
    tlevCounts = {}
    for i in range(0, len(taxonomies_list)):
        taxa = taxonomies_list[i]
        genus = taxa[5]
        if genus not in list(tlevCounts.keys()):
            tlevCounts[genus] = [1, i]
        else:
            tlevCounts[genus][0] += 1

    #find consensus (at the genus level?)
    maxValue = max(list(tlevCounts.values())[0])
    maxValKeys = [k for k, v in list(tlevCounts.items()) if v[0] == maxValue]

    if len(maxValKeys) == 1:
        consensus_index = tlevCounts[maxValKeys[0]][1]
        return taxonomies_list[consensus_index]
    else:
        # return get_lowest_common_taxonomy(opts, taxonomies_list)
        return get_lowest_common_taxonomy2(opts, taxonomies_list)
#

def assigntax2all(opts, asv, taxonomy_dict):

    #assign taxonomy to all ASVs, even host/contaminants/PhiX that we would normally remove
    best_hit_bitscore = "NA"
    best_hit_pident = "NA"
    best_hit_qcov = "NA"
    taxonomies_list = []
    # isPhiX = False
    # if asv.id == "Asv7":
        # pdb.set_trace() ###pdb

    #assign taxonomy to the "Keepers", which have at least one K or E hits (or N)
    if asv.designation == "K":
        if asv.get_best_hit_taxIDs("K") is not None:
            best_hit_bitscore = asv.get_best_hit_bitscore("K")
            best_hit_pident = asv.get_best_hit_pident("K")
            best_hit_qcov = asv.get_best_hit_qcov("K")
            #make a list of top taxonomies
            for txid in asv.get_best_hit_taxIDs("K"):
                if txid in taxonomy_dict:
                    taxonomies_list.append(taxonomy_dict[txid].split(";"))
                else:
                    print("K txid["+str(txid)+"] NOT in taxonomy_dict!", file=sys.stderr)

        elif asv.get_best_hit_taxIDs("E") is not None:
            best_hit_bitscore = asv.get_best_hit_bitscore("E")
            best_hit_pident = asv.get_best_hit_pident("E")
            best_hit_qcov = asv.get_best_hit_qcov("E")
            #make a list of top taxonomies
            #taxonomies_list = []
            for txid in asv.get_best_hit_taxIDs("E"):
                if txid in taxonomy_dict:
                    taxonomies_list.append(taxonomy_dict[txid].split(";"))
                else:
                    print("E txid["+str(txid)+"] NOT in taxonomy_dict!", file=sys.stderr)

    #assign taxonomy to the "Rejected" ASVs
    elif asv.designation == "R":
        if asv.get_best_hit_taxIDs("R") is not None:
            best_hit_bitscore = asv.get_best_hit_bitscore("R")
            best_hit_pident = asv.get_best_hit_pident("R")
            best_hit_qcov = asv.get_best_hit_qcov("R")
            # #deal with PhiX
            # if '10847' in asv.get_best_hit_taxIDs("R"):
                # #then we must assume this IS PhiX (other non-virus hits probably have PhiX contamination)
                # #because cramming other Kingdoms together results in "Unassigned"
                # #TODO: investigate/fix cases where best hits are also to bacteria and other viruses besides PhiX (these currently get "Unassigned")
                # taxonomies_list.append(taxonomy_dict['10847'].split(";"))
            # else:
                # #make a list of top taxonomies
                # for txid in asv.get_best_hit_taxIDs("R"):
                    # if txid in taxonomy_dict:
                        # taxonomies_list.append(taxonomy_dict[txid].split(";"))
                    # else:
                        # print("R txid["+str(txid)+"] NOT in taxonomy_dict!", file=sys.stderr)
            #make a list of top taxonomies
            for txid in asv.get_best_hit_taxIDs("R"):
                if txid in taxonomy_dict:
                    taxonomies_list.append(taxonomy_dict[txid].split(";"))
                else:
                    print("R txid["+str(txid)+"] NOT in taxonomy_dict!", file=sys.stderr)

        #only assign from "Undesignated" when it's the only choice
        elif "U" in asv.predesignation:
            best_hit_bitscore = asv.get_best_hit_bitscore("U")
            best_hit_pident = asv.get_best_hit_pident("U")
            best_hit_qcov = asv.get_best_hit_qcov("U")
            #make a list of top taxonomies
            for txid in asv.get_best_hit_taxIDs("U"):
                if txid in taxonomy_dict:
                    taxonomies_list.append(taxonomy_dict[txid].split(";"))
                else:
                    print("U txid["+str(txid)+"] NOT in taxonomy_dict!", file=sys.stderr)

    #lowest_common_taxa = get_lowest_common_taxonomy(opts, taxonomies_list, asv)
    lowest_common_taxa = get_lowest_common_taxonomy2(opts, taxonomies_list, asv)
    #adjusted_common_taxa = lowest_common_taxa

    #if all else fails... eg: 'N[(||NoTaxID|) ]'
    if len(lowest_common_taxa) == 0:
        lowest_common_taxa = ['k__Unassigned']

    taxonomy = asv.id+"\t"+';'.join(lowest_common_taxa)+"\t"+str(best_hit_bitscore)+"\t"+str(best_hit_pident)+"\t"+str(best_hit_qcov)
    if opts["--add_sizes"]:
        taxonomy = taxonomy+"\t"+str(asv.size)

    # if isPhiX:
        # taxonomy = None

    return taxonomy
#

def assigntax2keepers(opts, asv, taxonomy_dict):
    print("Warning: do you really want to use assigntax2keepers?")
    sys.exit(0)
    #we'll assign taxonomy to the "Keepers" only, which have at least one K or E hits (or N)
    best_hit_bitscore = "NA"
    taxonomies_list = []
    if asv.designation == "K":
        lowest_common_taxa = []
        if asv.get_best_hit_taxIDs("K") is not None:
            best_hit_bitscore = asv.get_best_hit_bitscore("K")
            #make a list of top taxonomies
            for txid in asv.get_best_hit_taxIDs("K"):
                taxonomies_list.append(taxonomy_dict[txid].split(";"))
            # lowest_common_taxa = get_consensus_taxonomy(opts, taxonomies_list)
            # lowest_common_taxa = get_lowest_common_taxonomy(opts, taxonomies_list, asv)
            lowest_common_taxa = get_lowest_common_taxonomy2(opts, taxonomies_list, asv)

        elif asv.get_best_hit_taxIDs("E") is not None: # and len(lowest_common_taxa) == 0:
            best_hit_bitscore = asv.get_best_hit_bitscore("E")
            #make a list of top taxonomies
            #taxonomies_list = []
            for txid in asv.get_best_hit_taxIDs("E"):
                taxonomies_list.append(taxonomy_dict[txid].split(";"))
            # lowest_common_taxa = get_consensus_taxonomy(opts, taxonomies_list)
            # lowest_common_taxa = get_lowest_common_taxonomy(opts, taxonomies_list, asv)
            lowest_common_taxa = get_lowest_common_taxonomy2(opts, taxonomies_list, asv)

        adjusted_common_taxa = lowest_common_taxa ## needed only because we commented out: # if len(lowest_common_taxa) > 0:

        #if all else fails... eg: 'N[(||NoTaxID|) ]'
        if len(adjusted_common_taxa) == 0:
            adjusted_common_taxa = ['Unassigned']

        taxonomy = asv.id+"\t"+';'.join(adjusted_common_taxa)+"\t"+str(best_hit_bitscore)
        if opts["--add_sizes"]:
            taxonomy = taxonomy+"\t"+str(asv.size)

    else:
        taxonomy = None

    return taxonomy
#

def assign_taxonomy(opts, asv_list, taxonomy_dict):
    '''
#Assignment Rules:
1 k <=86.1
2 p <=90 *interpolated
3 c <=100
4 o <140 *interpolated
5 f <180
6 g 268-180
    '''
    ### TODO: define these elsewhere and pass them in?
    # assignRules = {1:86.1, 2:90, 3:100, 4:140, 5:180, 6:268}
    ###

    assigned_taxonomy = []
    for asv in asv_list:
        # asvid = asv.id
        # if asvid == "Asv299": ###debuging###
            # pass
        if opts["--assign_all"]:
            #we'll assign taxonomy to everything - "Keepers" or "Rejected"
            taxonomy = assigntax2all(opts, asv, taxonomy_dict)
        else:
            #we'll assign taxonomy to the "Keepers" only, which have at least one K or E hits (or optionally N)
            taxonomy = assigntax2keepers(opts, asv, taxonomy_dict)

        if taxonomy is not None:
            assigned_taxonomy.append([taxonomy])

    return assigned_taxonomy
#

def main(args):
    opts = get_and_check_opts(args)
    opts["debug"] = False
    opts["--del_ssp"] = True
    opts['email'] = opts['-m'] # Always tell NCBI who you are

    ###
    if 1:
        if opts["debug"]:
            print(":parse_ASVs2filter_log", file=sys.stderr)

        asv_list = parse_ASVs2filter_log(opts)

        if opts["debug"]:
            print(":asv_list.sort", file=sys.stderr)
        asv_list.sort(key=lambda x: x.size, reverse=True)

        ### debug: save asv_list
        if 0:
            with open('asv_list.pickle', 'w') as f:
                print(":pickle.dump(asv_list.pickle)", file=sys.stderr)
                pickle.dump(asv_list, f)
                sys.exit(0)
        ###

    else:
        with open('asv_list.pickle') as f:
            print(":pickle.load(asv_list.pickle)", file=sys.stderr)
            asv_list = pickle.load(f)
    ###

    if opts["debug"]:
        print(":get_taxonomy(opts, asv_list)", file=sys.stderr)

    #get taxIDs from a local taxonomy DB and any newly-downloaded ones
    taxonomy_dict = get_taxonomy(opts, asv_list)

    # devel__print_blast_hit_summaries(opts, asv_list, taxonomy_dict)

    if opts["debug"]:
        print(":assign_taxonomy(opts, asv_list, taxonomy_dict)", file=sys.stderr)

    assigned_taxonomy = assign_taxonomy(opts, asv_list, taxonomy_dict)

    if opts["--add_sizes"]:
        header = "#ASVID\ttaxonomy\tbitscore\tper_id\tper_qcov\tsize"
    else:
        header = "#ASVID\ttaxonomy\tbitscore\tper_id\tper_qcov"

    if opts["-o"]:
        with open(opts["-o"], 'w') as outfile:
            outfile.write(header+"\n")
            for assndtax in assigned_taxonomy:
                outfile.write(assndtax[0]+"\n")
    else:
        print(header)
        for assndtax in assigned_taxonomy:
            print(assndtax[0])
#

if __name__ == "__main__":
    main(__doc__)