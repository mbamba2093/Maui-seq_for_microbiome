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

