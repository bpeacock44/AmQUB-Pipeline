#!/usr/bin/env python3
 
"""
assign_LCA_via_blast.py
 
Loads a BLAST output summary file created by [blast_top_hit_parser.py] and assigns taxonomy to the ASVs/query sequences.
 
Usage:
    assign_LCA_via_blast.py -i <file> -m <email@email.com> -o <file>
    assign_LCA_via_blast.py -h
 
Arguments:
    -i <file>     Path to the BLAST summary file (tophits_sum.txt)
    -m <email>    User's email for NCBI reporting
    -o <file>     Output file
 
Related scripts:
1) blast_top_hit_parser.py (creates the [tophits_sum.txt] file)
"""
 
from collections import defaultdict
from urllib.error import HTTPError
from docopt import docopt
from Bio import Entrez
import json
import time
import sys
import re
import os
import http.client
import argparse
import xmltodict
 
import re

class ASV:
    def __init__(self, asv):
        self.id = str(asv)
        self.hits = []
        self.best_hit_txIDS = []
        self.best_hit_bitscore = 0.0
        self.best_hit_pident = 0.0
        self.best_hit_qcov = 0.0
 
    def add_blast_hits(self, hitList):
        self.hits = hitList
        self.set_best_hits(hitList)
 
    def set_best_hits(self, hitList):
        if not hitList:
            return
        best_bitscore = 0.0
        for hit in hitList:
            parts = re.split(r'\|', hit)
            if len(parts) != 5:
                continue
            bitscore, pident, taxID, counts, qcov = parts
            bitscore = float(bitscore)
            if bitscore > best_bitscore:
                best_bitscore = bitscore
                self.best_hit_pident = float(pident)
                self.best_hit_qcov = float(qcov)
                self.best_hit_txIDS = [taxID]
            elif bitscore == best_bitscore:
                self.best_hit_txIDS.append(taxID)
        self.best_hit_bitscore = best_bitscore

def get_and_check_opts(args):
    parser = argparse.ArgumentParser(
        description="Loads a BLAST output summary file created by [blast_top_hit_parser.py] and assigns taxonomy to the ASVs/query sequences.",
        usage="%(prog)s -i <input_file> -m <email@email.com> -o <output_file>"
    )
    
    # Define expected arguments
    parser.add_argument('-i', type=str, required=True, help='Filepath to the BLAST summary file (tophits_sum.txt)')
    parser.add_argument('-m', type=str, required=True, help='User\'s email for NCBI reporting')
    parser.add_argument('-o', type=str, required=True, help='Output file')
 
    # Parse arguments
    opts = parser.parse_args(args)
    
    optionsFail = False
 
    # Validate the email format
    if not re.match(r"[^@]+@[^@]+\.[^@]+", opts.m):
        print("Invalid email format")
        optionsFail = True
    
    # Validate file paths
    if not os.path.isfile(opts.i):
        print(f"Input file '{opts.i}' does not exist")
        optionsFail = True
 
    # Exit if any validation fails
    if optionsFail:
        sys.exit(1)  # Exit with a non-zero status to indicate an error
    
    return vars(opts)  # Convert Namespace to dictionary
 
def parse_asv_file(file_path):
    asvs = {}
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()
                if not line or line.isspace():
                    continue
                try:
                    asv_id, hit_data_str = line.split('\t')
                except ValueError:
                    print(f"Skipping line due to incorrect format: {line}")
                    continue
                # No need to split the ASV ID or handle size
                asv_name = asv_id
                # Process the hit data string
                hit_data_str = hit_data_str.strip('[]')
                if hit_data_str == '':
                    hits = []
                else:
                    hits = hit_data_str.split(') (')
                    hits = [hit.replace('(', '').replace(')', '') for hit in hits]
                # Create an ASV instance and add the hits
                asv = ASV(asv_name)
                asv.add_blast_hits(hits)
                asvs[asv_name] = asv
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
    return asvs

 
def parse_xml_file(taxa_xml_file):
    """
    Parse the XML file and return a dictionary of taxon IDs to taxonomic lineages.
    @param taxa_xml_file: Path to the XML file
    @return: dict[taxonIDs]=taxonomy
    """
    d = {}
    taxa_level_prefix = "k__ p__ c__ o__ f__ g__ s__ u__".split()
    ranks = "superkingdom phylum class order family genus species".split()  # Removed subspecies
    if not os.path.isfile(taxa_xml_file):
        print(f":  ** File not found: [{taxa_xml_file}] **")
        return d    
    try:
        with open(taxa_xml_file) as fh:
            rootD = xmltodict.parse(fh.read())
    except Exception as e:
        print(f":  ** Error reading file [{taxa_xml_file}]: {e} **")
        return d
    if isinstance(rootD, dict) and isinstance(rootD.get('TaxaSet', {}), dict):
        taxon_objectsL = rootD['TaxaSet'].get('Taxon', [])
    else:
        print(f":  ** Nothing to parse in [{taxa_xml_file}] **")
        return d    
    if isinstance(taxon_objectsL, dict):
        taxon_objectsL = [taxon_objectsL]    
    for TaxonD in taxon_objectsL:
        TaxId = TaxonD.get('TaxId', 'Unknown')
        ScientificName = TaxonD.get('ScientificName', 'Unknown')
        taxonL = TaxonD.get('LineageEx', {}).get('Taxon', [])
        if isinstance(taxonL, dict):
            taxonL = [taxonL]        
        taxonomyD = {rank: None for rank in ranks}
        genus, species = "", ""
        species_found = False        
        # Initialize taxonomy list
        taxonomy = []
        last_defined_taxa = ""        
        for taxlevD in taxonL:
            rank = taxlevD.get('Rank', 'Unknown')
            if rank in taxonomyD:
                taxonomyD[rank] = taxlevD.get('ScientificName', 'Unknown')        
        for i, rank in enumerate(ranks):
            if rank == 'genus':
                if taxonomyD[rank] is not None:
                    genus = taxonomyD[rank]
                    species_name = re.sub(f"{genus} ", "", ScientificName).split(" ")[0]
                    species = f"{genus}_{species_name}"
                    taxonomyD['species2'] = species
                else:
                    taxonomyD['species2'] = None
            elif rank == 'species':
                if taxonomyD[rank] is not None:
                    species_found = True
                    species = taxonomyD[rank]
                    taxonomyD['species'] = species
                else:
                    taxonomyD['species'] = taxonomyD.get('species2', None)
            # Skip subspecies
            if taxonomyD.get(rank) is not None:
                taxonomy_rank = re.sub(" ", "_", taxonomyD[rank])
                taxonomy.append(taxa_level_prefix[i] + taxonomy_rank)
                last_defined_taxa = taxonomy_rank
            else:
                taxonomy.append(taxa_level_prefix[i] + "unclassified_" + last_defined_taxa)        
        lineage = ";".join(taxonomy).replace(" ", "").replace("=", ".").replace("'", "")
        d[TaxId] = lineage        
        AkaTaxIds = TaxonD.get('AkaTaxIds', {}).get('TaxId', '')
        if AkaTaxIds:
            AkaTaxIds_list = re.split(",", AkaTaxIds)
            for AkaTaxId in AkaTaxIds_list:
                d[AkaTaxId] = lineage    
    return d
 
def epost_taxonIDs(opts, new_taxonIDs_list):
    #TODO: define these earlier? somewhere else?
    #set some options
    opts['taxa_xml_file'] = "taxa_01." # TODO: define directory as /tmp
    opts['retmax'] = 1000
    opts['retmode'] = "xml"
    opts['max_attempts'] = 3
    opts['count'] = len(new_taxonIDs_list)
 
    if opts['count'] > 0:
        print(":  ",str(opts['count'])+" taxIDs will be ePost'ed to NCBI")
        #epost taxon ID list
        opts['db'] = "taxonomy";
        Entrez.email = opts['email']
        epost_results = Entrez.read(Entrez.epost(opts['db'], id=",".join(new_taxonIDs_list)))
        #get WebEnv and QueryKey from epost_results
        opts['WebEnv']   = epost_results["WebEnv"]
        opts['QueryKey'] = epost_results["QueryKey"]
 
def get_next_xml_idx(opts):
    """
    Get the next available index for XML files.
    @param opts: Options dictionary
    @return: Next XML file index
    """
    idx = 1
    while os.path.isfile(opts['taxa_xml_file'] + str(idx) + ".xml"):
        idx += 1
    return idx 
 
def download_eposted_taxonIDs_to_XML(opts):
    """
    Download taxonomies in XML format and save them to files.
 
    @param opts: Options dictionary containing configuration and parameters
    """
    opts['xml_files'] = []  # Initialize the list to store names of XML files
    idx = get_next_xml_idx(opts)  # Get the index for naming XML files
    
    for start in range(0, opts['count'], opts['retmax']):
        out_file = f"{opts['taxa_xml_file']}{idx}.xml"
        with open(out_file, 'w') as out_handle:
            end = min(opts['count'], start + opts['retmax'])
            num2fetch = str(end - start)
            attempt = 1
            
            while attempt <= opts['max_attempts']:
                print(f": Efetching taxa. Attempt {attempt}")
                try:
                    fetch_handle = Entrez.efetch(
                        db="taxonomy",
                        rettype="xml",
                        retmode=opts['retmode'],
                        retstart=start,
                        retmax=opts['retmax'],
                        webenv=opts['WebEnv'],
                        query_key=opts['QueryKey']
                    )
                    
                    data = b''
                    while True:
                        chunk = fetch_handle.read(1024)
                        if not chunk:
                            break
                        data += chunk
                    
                    fetch_handle.close()
                    out_handle.write(data.decode('utf-8'))
                    opts['xml_files'].append(out_file)
                    idx += 1
                    print(f": Taxa fetched successfully for {num2fetch} taxa.")
                    break
                
                except (HTTPError, http.client.IncompleteRead) as err:
                    print(f": Received error from server {err}")
                    attempt += 1
                    time.sleep(15 * attempt)
                
                except Exception as e:
                    print(f": Error occurred: {e}")
                    raise
 
    print(f": XML files downloaded: {', '.join(opts['xml_files'])}")
 
def get_taxonomies_from_XML_files(opts):
    """
    Retrieve taxonomies from downloaded XML files.
 
    @param opts: Options dictionary containing 'xml_files' key with list of XML file paths
    @return: dict[taxonIDs]=taxonomic_lineages
    """
    all_new_taxonomies_dict = {}
    
    for out_file in opts['xml_files']:
        print(f": Parsing XML file [{out_file}]")
        try:
            # Parse the XML file and update the taxonomy dictionary
            new_taxIDs_dict = parse_xml_file(out_file)
            all_new_taxonomies_dict.update(new_taxIDs_dict)
        except Exception as e:
            print(f"** Error parsing file [{out_file}]: {e} **")
    
    return all_new_taxonomies_dict
 
def get_taxonIDs_list(asv_list):
    """
    Retrieve a list of unique taxon IDs from a list of ASV objects.
 
    @param asv_list: List of ASV objects
    @return: List of unique taxon IDs
    """
    taxIDs = set()  # Use a set to automatically handle duplicates
    
    for asv in asv_list:
        taxIDs.update(asv.best_hit_txIDS)  # Access attribute directly
    
    return list(taxIDs)  # Convert set back to list
 
def update_taxonomy(self, new_taxonomies_dict):
    """
    Update the ASV object with new taxonomy information.
    @param new_taxonomies_dict: Dictionary of taxonomies
    """
    for taxID in self.best_hit_txIDs():
        if taxID in new_taxonomies_dict:
            self.best_hit_taxonomy = new_taxonomies_dict[taxID] 
 
def get_taxonomy(opts, asv_list):
    """
    Retrieve and update taxonomic information based on ASV data.
 
    @param opts: Options dictionary containing configuration and file paths
    @param asv_list: Dictionary of ASV objects returned by parse_asv_file
    @return: Dictionary of taxonomic information indexed by taxon IDs
    """
    # Extract taxon IDs from the ASV list
    new_taxonIDs_list = get_taxonIDs_list(asv_list.values())
    
    new_taxonomies_dict = {}
    
    if new_taxonIDs_list:
        # Post taxon IDs to NCBI and download the resulting XML files
        epost_taxonIDs(opts, new_taxonIDs_list)
        download_eposted_taxonIDs_to_XML(opts)
        
        # Retrieve taxonomies from downloaded XML files
        new_taxonomies_dict = get_taxonomies_from_XML_files(opts)
    
    return new_taxonomies_dict
 
def get_lowest_common_taxonomy(opts, taxonomies_list, asv):
    """
    Retrieve the lowest common taxonomy from a list of taxonomic lineages.
 
    @param opts: Options dictionary
    @param taxonomies_list: List of taxonomic lineages (lists of taxonomic levels)
    @param asv: ASV object (not used directly in function but included in parameters for compatibility)
    @return: List of lowest common taxonomy levels
    """
    if not taxonomies_list:
        return []
 
    # Get the common taxonomy levels
    lowest_common_taxa = taxonomies_list[0]
    for taxa in taxonomies_list[1:]:
        lowest_common_taxa = [
            level for level, taxon in zip(lowest_common_taxa, taxa)
            if level == taxon and "__unclassified" not in level
        ]
 
    # If no common taxonomy found, find the most common taxonomy levels
    if not lowest_common_taxa:
        dc = defaultdict(lambda: defaultdict(int))
        minrange = min(map(len, taxonomies_list))
        
        # Count occurrences of cumulative taxonomies
        for level in range(1, minrange + 1):
            for txa in taxonomies_list:
                cumTaxa = ';'.join(txa[:level])
                dc[level][cumTaxa] += 1
 
        # Find the most common taxonomy at each level
        for level in range(1, minrange + 1):
            max_taxa = max(dc[level], key=dc[level].get)
            max_taxa_count = dc[level][max_taxa]
            if sum(count == max_taxa_count for count in dc[level].values()) >= 1:
                if level == 1:
                    all_taxa_at_max_count = sorted(dc[level].keys())
                    lowest_common_taxa = ["_OR_".join(all_taxa_at_max_count)]
                else:
                    back_level = level - 1
                    max_taxa_at_back_level = max(dc[back_level], key=dc[back_level].get)
                    lowest_common_taxa = max_taxa_at_back_level.split(";")
                break
 
    return lowest_common_taxa
 
def assign_taxonomy(opts, asv_list, taxonomy_dict):
    assigned_taxonomy = []
    for asv in asv_list.values():
        taxonomies_list = []
        # Get best hit details
        best_hit_taxIDs = asv.best_hit_txIDS
        if best_hit_taxIDs:
            best_hit_bitscore = asv.best_hit_bitscore
            best_hit_pident = asv.best_hit_pident
            best_hit_qcov = asv.best_hit_qcov
            # Make a list of top taxonomies
            for txid in best_hit_taxIDs:
                if txid in taxonomy_dict:
                    taxonomies_list.append(taxonomy_dict[txid].split(";"))
                else:
                    print(f"txid[{txid}] NOT in taxonomy_dict!")
            lowest_common_taxa = get_lowest_common_taxonomy(opts, taxonomies_list, asv)
        # If no common taxonomy found, use 'Unassigned'
        if not lowest_common_taxa:
            lowest_common_taxa = ['k__Unassigned']
        taxonomy = f"{asv.id}\t{';'.join(lowest_common_taxa)}\t{best_hit_bitscore}\t{best_hit_pident}\t{best_hit_qcov}"
        assigned_taxonomy.append([taxonomy])
    return assigned_taxonomy
 
def main(args):
    opts = get_and_check_opts(args)

    # Always tell NCBI who you are
    opts['email'] = opts['m']
    
    # Parse the ASV file
    asv_list = parse_asv_file(opts['i'])
    
    # Get taxIDs from a local taxonomy DB and any newly-downloaded ones
    taxonomy_dict = get_taxonomy(opts, asv_list)
    
    # Assign taxonomy based on the retrieved taxIDs
    assigned_taxonomy = assign_taxonomy(opts, asv_list, taxonomy_dict)
    
    # Define the header for the output file
    header = "#ID\ttaxonomy\tbitscore\tper_id\tper_qcov"
    
    # Write the results to the output file
    with open(opts["o"], 'w') as outfile:
        outfile.write(header + "\n")
        for assndtax in assigned_taxonomy:
            outfile.write("\t".join(assndtax) + "\n")
 
if __name__ == "__main__":
    main(sys.argv[1:])
