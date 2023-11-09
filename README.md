# Maui-seq_for_microbiome

This script is modification of the method described in MAUI-seq (Moeskjar et al. 2020 Mol Ecol Resour https://doi.org/10.1111/1755-0998.13294), aiming to simplify the process and enable its application to large-scale microbial community data. In this analysis, the amplified product profile in the sample is evaluated by counting the number of random sequence tags (unique molecular identifiers, UMIs) incorporated into the primers that are observed in a particular amplified product.

When using this script, you will need to set the paths for PEAR and BLAST+ and prepare the BLAST reference database along with the corresponding taxonomic information. The required files used in Bamba et al. 2023 are stored in the "test_data" directory. Please refer to that directory for the necessary files.

## MAUI_main.py
This script serves as a simple executable file that imports the "MAUI_modules.py" module.

The usage of this script is as follows:

```python MAUI_main.py MAUI_config.ini```

## MAUI_modules.py
Functions stored in this file are imported by the main file.

The code is a Python script for executing the MAUI-seq data analysis pipeline. 

The pipeline automates a series of steps including preprocessing, assembly, counting, ID conversion, creation of a taxonomic database through BLAST, and assignment of taxonomy based on database IDs.



## MAUI_config.ini
This configuration file provides parameters for PEAR and BLAST and the paths to the various data.

[Main_configure]

data_dir = /path/to/data_directory/ # PATH to directory containing .fastq files.

out_dir = /path/to/output_directory/ # PATH to output directory.

prefix = First_analysis # Prefixes for result files.

need_blast = True # if False, BLAST search will not be performed.

taxa_db = /path/to/taxa_db/ # PATH to output taxa database. If not found, it is generated.

seq_db = /path/to/seq_db/ # PATH to output taxa database. If not found, it is generated.

[Maui_configure]

marker_gene = V5 # Please see the MAUI_modules.Config.maui_gene(). This information is used for detect primer length. 

UMI_len = 12 # How many bases are in the UMI sequence upstream of the primer sequence.

min_ave_reads_number = 0.1 # Minimum number of UMI sequences allowed for all sample averages.


[Pear_configure]

min_overlap = 10

min_assembly_length = 300

min_trim_length = 200

quality_threashold = 20

number_of_threads = 42

[Blast_configure]

blast_db = /path/to/blast_db # PATH to BLAST+ database (output makeblastdb).

blast_threads = 42

max_target_seqs = 5

taxonomic_database = /path/to/reference_taxa_db # PATH to taxonomic information of BLAST+ database,

