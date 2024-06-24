# Qiime2 Implementation
If you would like to assign taxonomy using a curated database instead of using the NCBI nt database, Qiime2 provides excellent tools for doing so. 

Briefly, Qiime2 recommends using a classifier to assign taxonomy. Ideally, this classifier will be trained specifically for your particular data and using the database of your choice. (See Qiime2 website for details and tools to do this.) 

As of 3/5/24, Qiime2 also provides four pre-trained classifiers - one trained using the SILVA database, one trained using SILVA (just 515F/806R region of sequences), one trained using the Greengenes database, and one trained using Greengenes (just 515F/806R region of sequences). These are very easy to use.

In our intial tests of universal assay data, the SILVA classifier worked fairly well to provide a similar prediction as a manual BLAST might lead one to provide. However, we found that in many cases it predicted taxonomy down to genus/species level even though a BLAST search of the same sequence indicated that the taxonomy couldn't be resolved past family level. This may or may not be an issue, depending on your particular analysis. (Note that Silva does not curate to Species level and at some point Qiime will no longer proivde species-level predictions as a result.)

## Using Qiime2 for Taxonomic Assignment of ASVs
If you complete this pipeline through mbio_part3.sh, then you will have an initial ASV table in biom format and a fasta file of your ASVs:

- path/to/out_directory/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta
- path/to/out_directory/asvs/asv_table_01.biom

You can use these, for example, as follows within a conda environment created with Qiime2:

```sh
# Create conda environment with Qiime2:
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-osx-conda.yml
conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-osx-conda.yml
rm qiime2-amplicon-2024.2-py38-osx-conda.yml
conda activate qiime2-amplicon-2024.2

# Download SILVA classifier
wget https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza

# Create a qza file of your ASV sequences
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path <output directory>/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta \
  --output-path sequences.qza

# Run classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads sequences.qza \
  --o-classification taxonomy.qza

# View results
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

Further analyses can be completed using the .biom file and this taxonomy.qzv file using Qiime2 tools.
