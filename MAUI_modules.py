#!/usr/bin/env python3
import sys
import os
import sys
import glob
import numpy as np
import pandas as pd
import subprocess as sp
import configparser
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import time
import itertools


class Module_check:
    #This pipeline for MAUI-seq analysis utilizes BLAST+ and PEAR in the process of analysis.
    #Please ensure that the necessary software is available for use.

    def __init__(self):
        def modules():
            blast = shutil.which("blastn")
            pear = shutil.which("pear")
            if blast == None:
                print("ERROR: BLAST+ is not found. Please make path to blastn function as blastn.")
                sys.exit()

            if pear == None:
                print("ERROR: PEAR is not found. Please make path to PEAR as pear.")
                sys.exit()

        modules()
        print("Modules were ready.")





class Config:
    #In this analysis, we use an .ini file as the configuration file.
    #Please create a file with the necessary items filled in.

    config_file = sys.argv[1]
    config = configparser.ConfigParser()
    config.read(config_file)

    def data_dir(self):
        return self.config["Main_configure"].get("data_dir")

    def out_dir(self):
        return self.config["Main_configure"].get("out_dir")

    def prefix(self):
        return self.config["Main_configure"].get("prefix")

    def need_blast(self):
        return self.config["Main_configure"].get("need_blast")

    def taxa_db(self):
        return self.config["Main_configure"].get("taxa_db")

    def seq_db(self):
        return self.config["Main_configure"].get("seq_db")

    def maui_gene(self):
        marker = self.config["Maui_configure"].get("marker_gene")
        marker_list = {
        "V5": {"f_primer_len":19,"r_primer_len": 19},
        "V3": {"f_primer_len":17,"r_primer_len": 21},
        "nosZC1": {"f_primer_len":17,"r_primer_len": 19},
        "nosZC2": {"f_primer_len":22,"r_primer_len": 20},
        "rpoB": {"f_primer_len":18,"r_primer_len": 17},
        "AMF": {"f_primer_len":23,"r_primer_len": 21}}
        if not marker in marker_list.keys():
            print("Error: Marker gene is not supported in this script.\nPlease provide primer length to config_check.maui_gene.")
            sys.exit()
        else:
            return marker_list[marker]

    def maui_params(self):
        UMI_len = int(self.config["Maui_configure"].get("UMI_len"))
        threashold = float(self.config["Maui_configure"].get("min_ave_reads_number"))
        params = {"UMI_len": UMI_len, "threashold": threashold}
        return params

    def pear_params(self):
        min_overlap = self.config["Pear_configure"].get("min_overlap")
        min_assembly_length = self.config["Pear_configure"].get("min_assembly_length")
        min_trim_length = self.config["Pear_configure"].get("min_trim_length")
        quality_threashold = self.config["Pear_configure"].get("quality_threashold")
        number_of_threads = self.config["Pear_configure"].get("number_of_threads")
        params = f" -v {min_overlap} -n {min_assembly_length} -t {min_trim_length} -q {quality_threshold} -j {number_of_threads}"
        return params

    def blast_params(self):
        blast_db = self.config["Blast_configure"].get("blast_db")
        blast_threads = self.config["Blast_configure"].get("blast_threads")
        max_target_seqs = self.config["Blast_configure"].get("max_target_seqs")
        params = f" -db {blast_db} -num_threads {blast_threads} -max_target_seqs {max_target_seqs} -outfmt 10 -task blastn"
        return params

    def taxonomic_database(self):
        return self.config["Blast_configure"].get("taxonomic_database")





class PreProcess:
    def __init__(self, config):
        # Create the output directory specified in the config
        os.makedirs(config.out_dir(), exist_ok = True)

        # Check if Blast is needed according to the config
        if config.need_blast():
            # Check if the sequence database file exists
            neworupdate = os.path.isfile(config.seq_db())

            # If the sequence database file doesn't exist, create empty databases
            if neworupdate != True:
                seq_DB = dict()
                with open(config.seq_db(), "wb") as f:
                    pickle.dump(seq_DB, f)
                taxa_DB = dict()
                with open(config.taxa_db(), "wb") as f:
                    pickle.dump(taxa_DB, f)

            # If the sequence database file exists, create backups and load existing databases.
            else:
                copytime = time.time()
                if os.path.isfile(config.taxa_db()) == False:
                    taxa_DB = dict()
                    with open(config.taxa_db(), "wb") as f:
                        pickle.dump(taxa_DB, f)

                # Load the existing taxa database.
                with open(config.taxa_db(), "rb") as f:
                    taxa_DB = pickle.load(f)

                # Create a backup of the taxa database with a timestamp.
                copiedtaxa_db = config.taxa_db().split(".pkl")[0] + str(round(copytime)) + ".pkl"
                with open(copiedtaxa_db, "wb") as f:
                    pickle.dump(taxa_DB, f)

                with open(config.seq_db(), "rb") as f:
                    seq_DB = pickle.load(f)

                copiedseq_db = config.seq_db().split(".pkl")[0] + str(round(copytime)) + ".pkl"
                with open(copiedseq_db, "wb") as f:
                    pickle.dump(seq_DB, f)

        print("PreProcess were finished.")





class Pear:
    def __init__(self, config):
        # Get the data and output directories from the configuration.
        data_dir = config.data_dir()
        out_dir = config.out_dir()

        # Retrieve the paths of the forward and reverse fastq files using glob and sort them.
        fastq_f = sorted(glob.glob(data_dir + "*_R1*.fastq.gz"))
        fastq_r = sorted(glob.glob(data_dir + "*_R2*.fastq.gz"))

        ### Fastq file check ###
        # Check if the numbers of forward and reverse fastq files are the same.
        if not len(fastq_f) == len(fastq_f):
            print("Error: Numbers of fastq files were different between R1 and R2.")
            sys.exit()

        check_list = []
        for i, j in zip(fastq_f, fastq_r):
            # Extract the sample name from the file paths.
            temp_i = i.split("/")[-1].split("_R1")[0]
            temp_j = j.split("/")[-1].split("_R2")[0]

            # Check if the sample names match between the forward and reverse fastq files.
            check_list.append(temp_i == temp_j)

        if not all(check_list):
            print("Error: Fastq files could not be paired")
            sys.exit()
        ### Fastq file check fin ###

        # Retrieve the PEAR parameters from the configuration.
        pearparams = config.pear_params()

        # Perform PEAR assembly for each pair of forward and reverse fastq files.
        for i, j in zip(fastq_f, fastq_r):
            # Extract the sample name from the file path.
            output_filename = i.split("/")[-1].split("_")[0]
            output_filename = out_dir + output_filename

            # Build the PEAR command and execute it using subprocess.
            Pear_call = f"pear -f {i} -r {j} -o {output_filename} {pearparams}"
            sp.call(Pear_call, shell = True)

        # Create directories for storing the PEARed and unPEARed data.
        os.makedirs(f"{out_dir}Data_MAUI_PEARed", exist_ok = True)
        os.makedirs(f"{out_dir}Data_MAUI_unPEARed", exist_ok = True)

        # Move the assembled fastq files to the "Data_MAUI_PEARed" directory.
        assembled_file_list = sorted(glob.glob(out_dir + "*assembled.fastq"))
        index_list = [i.split("/")[-1].split(".assembled.fastq")[0] for i in assembled_file_list]
        assembled_read_count = []
        for i in assembled_file_list:
            # Count the number of reads in each assembled fastq file.
            assembled_read_count.append(sum([1 for _ in open(i)]))
            shutil.move(i, f"{out_dir}Data_MAUI_PEARed")

        # Move the unPEARed fastq files to the "Data_MAUI_unPEARed" directory.
        unassembled_file_list = sorted(glob.glob(out_dir + "*.fastq"))
        for i in unassembled_file_list:
            shutil.move(i, f"{out_dir}Data_MAUI_unPEARed")

        # Count the number of reads in each discarded fastq file.
        unassembled_read_count = []
        discarded_file_list = sorted(glob.glob(out_dir + "Data_MAUI_unPEARed/*discarded.fastq"))
        for i in discarded_file_list:
            unassembled_read_count.append(sum([1 for _ in open(i)]))

        # Create a DataFrame to store the read numbers and save it as a CSV file.
        read_number_df = pd.DataFrame(index = index_list, columns = ["Assembled", "Unassembled"])
        read_number_df["Assembled"] = assembled_read_count
        read_number_df["Unassembled"] = unassembled_read_count
        read_number_df.to_csv(out_dir + "Read_number.csv")

        print("Pair-end read assembly with pear was finished.")





class MauiCounting:
    def __init__(self, config):
        ### Configure check ###
        # Get the output directory and prefix from the configuration.
        out_dir = config.out_dir() + "Data_MAUI_PEARed/"
        prefix = config.prefix()

        # Retrieve the lengths of the forward and reverse primers from the Maui gene configuration.
        f_primer_len = config.maui_gene()["f_primer_len"]
        r_primer_len = config.maui_gene()["r_primer_len"]

        # Get the list of fastq filenames in the output directory.
        fastq_filename_list = glob.glob(out_dir + "*.fastq")

        # Retrieve the UMI length and minimum average reads number from the Maui parameters configuration.
        UMI_len = config.maui_params()["UMI_len"]
        min_ave_reads_number = config.maui_params()["threashold"]
        ### Configure check fin ###

        # Create a temporary dictionary to store the counts of each unique sequence in each fastq file.
        dict_temp = dict()
        for j in fastq_filename_list:
            with open(j) as fastq_file:
                # Extract the amplified sequences and remove the UMI and primers from each sequence.
                seq_list_temp = [str(i.seq[(UMI_len + f_primer_len): (len(i.seq) - r_primer_len)]) for i in SeqIO.parse(fastq_file, "fastq")]

                # Extract unique amplified sequences.
                seq_list_temp = set(seq_list_temp)

                # Counting MAUI for each sample.
                dict_temp[j] = dict()
                for i in seq_list_temp:
                    dict_temp[j][i] = list()
                for i in SeqIO.parse(j, "fastq"):
                    UMI = str(i.seq[0:UMI_len]) # Extract UMI nucleotide sequence
                    dict_temp[j][str(i.seq[(UMI_len + f_primer_len): (len(i.seq) - r_primer_len)])].append(UMI) # Add UMI nucleotide sequence to amplified sequence
                for i in dict_temp[j].keys():
                    dict_temp[j][i] = len(set(dict_temp[j][i])) # Count UMIs for each amplified sequence

        # Create a list of all unique sequences from the dictionary.
        seq_list = list()
        for j in fastq_filename_list:
            seq_list.extend(list(dict_temp[j].keys()))
        seq_list = list(set(seq_list))

        # Calculate the sum of counts for each sequence across all fastq files.
        sum_list = list()
        for i in seq_list:
            sum_temp = 0
            for j in fastq_filename_list:
                if i in dict_temp[j].keys():
                    sum_temp = sum_temp + dict_temp[j][i]
            sum_list.append(sum_temp)

        # Filter out sequences that have a sum of counts below the threshold.
        seq_list_filt = [seq_list[i] for i in range(0, len(sum_list)) if sum_list[i] < len(fastq_filename_list)*min_ave_reads_number]

        # Remove the filtered sequences from the dictionary.
        for i in seq_list_filt:
            for j in fastq_filename_list:
                if i in dict_temp[j].keys():
                    dict_temp[j].pop(i)

        # Save the temporary dictionary as a pickle file.
        with open(f"{out_dir}{prefix}_temporary.pkl", "wb") as f:
            pickle.dump(dict_temp, f)

        # Create a list of DataFrames for each fastq file's counts.
        temp_list = list()
        for i in dict_temp.keys():
            temp_list.append(pd.DataFrame(pd.Series(dict_temp[i], name = i)).T)

        # Concatenate the DataFrames into a single DataFrame.
        comp_df = pd.concat(temp_list, sort = True)
        comp_df = comp_df.fillna(0)

        # Calculate the sum of counts for each column (sequence) and set the index.
        sum_list = [sum(comp_df.iloc[: ,i]) for i in range(0, comp_df.shape[1])]
        index_list = [i.split("/")[-1].split(".")[0] for i in comp_df.index]
        comp_df.index = index_list

        # Add a row for the total counts and sort the DataFrame by the total counts column.
        comp_df.loc["Total"] = sum_list
        comp_df.sort_values(by = "Total", axis = 1, ascending = False , inplace = True)

        # Assign unique IDs to the sequences and save them as a FASTA file.
        seq_id_num = 0
        seq_list = list()
        seq_id_list = list()
        for i in comp_df.columns:
            seq_id = "seq_" + str(seq_id_num)
            seq = Seq(i)
            seq_r = SeqRecord(seq, id = seq_id, description = "")
            seq_id_num += 1
            seq_list.append(seq_r)
            seq_id_list.append(seq_id)
        SeqIO.write(seq_list, f"{out_dir}{prefix}.fas", "fasta")

        # Update the column names in the DataFrame with the sequence IDs and save it as a CSV file.
        comp_df.columns = seq_id_list
        comp_df.to_csv(f"{out_dir}{prefix}.csv")

        # Create a "fastq" directory and move the assembled fastq files into it.
        os.makedirs(f"{out_dir}fastq", exist_ok = True)
        assembled_file_list = glob.glob(out_dir + "*assembled.fastq")
        for i in assembled_file_list:
            shutil.move(i, f"{out_dir}fastq")

        print("MAUI counting was finished.")



class SeqID2DBid:
    def __init__(self, config):
        # Get the output directory, prefix, and file names from the configuration.
        out_dir = config.out_dir() + "Data_MAUI_PEARed/"
        prefix = config.prefix()
        comp_df_name = f"{out_dir}{prefix}.csv"
        fasta_name = f"{out_dir}{prefix}.fas"
        seq_DB_name = config.seq_db()
        need_blast = config.need_blast()
        blast_params = config.blast_params()
        taxa_DB_name = config.taxa_db()

        def seqIDupdate():
            # Read the composition table CSV and check the IDs of the output fasta file and composition table.
            comp_df = pd.read_csv(comp_df_name, index_col = 0)
            k = len(list(SeqIO.parse(fasta_name, "fasta")))
            id_list = [None]*k
            k1 = 0
            for i in SeqIO.parse(fasta_name, "fasta"):
                id_list[k1] = i.id
                k1 += 1
            if not all(id_list == comp_df.columns):
                print("IDs were not matched between the output fasta file and composition table.")
                sys.exit()

            # Load the sequence database from the file.
            with open(seq_DB_name, "rb") as f:
                seq_DB = pickle.load(f)

            # Check if the sequence database is empty (new) or has existing entries (update).
            if len(seq_DB.keys()) == 0:# New seq_DB
                startDBid = 0
            else:# Update seq_DB
                startDBid = int(np.max(np.array([int(i.split("_")[1]) for i in seq_DB.values()])))

            non_ovlap_fasta = list()
            fasta_list = list()
            fasta_id_list = list()
            # Process each sequence in the output fasta file.
            for i in SeqIO.parse(fasta_name, "fasta"):
                if not str(i.seq) in seq_DB.keys():
                    startDBid += 1
                    temp_DBid = "DB_" + str(startDBid)
                    seq_DB[str(i.seq)] = temp_DBid
                    i.id = temp_DBid
                    i.description = ""
                    non_ovlap_fasta.append(i)
                    fasta_list.append(i)
                    fasta_id_list.append(i.id)
                    seq_DB[str(i.seq)] = i.id
                else:
                    temp_DBid = seq_DB[str(i.seq)]
                    i.id = temp_DBid
                    i.description = ""
                    fasta_list.append(i)
                    fasta_id_list.append(i.id)

            # Write the non-overlapping fasta sequences to a separate file for BLAST search.
            SeqIO.write(non_ovlap_fasta, f"{out_dir}{prefix}_non_overlap.fas", "fasta")
            # Write the renamed fasta sequences to a separate file.
            SeqIO.write(fasta_list, f"{out_dir}{prefix}_rename.fas", "fasta")

            # Update the column names in the composition table and save it as a CSV file.
            comp_df.columns = fasta_id_list
            comp_df = comp_df[:-1]
            comp_df.to_csv(f"{out_dir}{prefix}_rename.csv")

            # Save the updated sequence database to the file.
            with open(seq_DB_name, "wb") as f:
                pickle.dump(seq_DB, f)

        def blast():
            # Run BLAST with the non-overlapping fasta file as the query.
            Blast_query = f"{out_dir}{prefix}_non_overlap.fas"
            Blast_command = f"blastn -query {Blast_query} -out {out_dir}{prefix}_BLAST.csv {blast_params}"
            sp.call(Blast_command, shell = True)

        def maketaxadb():
            # Read the BLAST output file and check if it contains results.
            Blast_out_name = f"{out_dir}{prefix}_BLAST.csv"
            with open(Blast_out_name) as f:
                Blast_out_check = f.readlines()

            if not len(Blast_out_check) == 0:
                # Read the BLAST output and the taxonomic database.
                Blast_out = pd.read_csv(Blast_out_name, index_col = 0, header = None)
                taxonomic_database_name = config.taxonomic_database()
                with open(taxonomic_database_name, "rb") as f:
                    taxonomic_database = pickle.load(f)

                # Load the existing taxa database.
                with open(taxa_DB_name, "rb") as f:
                    taxa_DB = pickle.load(f)
                taxa_DB_filesize1 = sys.getsizeof(taxa_DB)
                DBid_list = set(Blast_out.index)

                # Process each DB ID from the BLAST output.
                for i in DBid_list:
                    _temp_taxa = Blast_out.loc[Blast_out.index == i, ]
                    if len(_temp_taxa) == 0:
                        _refseq = ("not_assigned")
                        _taxa = "not_assigned"
                    else:
                        _temp_taxa_max = list(_temp_taxa[1][_temp_taxa[2] == max(_temp_taxa[2])])
                        _temp_boolean = [taxonomic_database[j]["genus"] != "NotaAssigned" for j in _temp_taxa_max]
                        _temp_taxa_filtered = list(itertools.compress(_temp_taxa_max, _temp_boolean))
                        if not _temp_taxa_filtered:
                            _refseq = _temp_taxa_max[0]
                        else:
                            _refseq = _temp_taxa_filtered[0]
                        if _refseq in taxonomic_database.keys():
                            _taxa = taxonomic_database[_refseq]
                        else:
                            assert False, ("check the blast results and database, there are non-matched sequences")
                    taxa_DB[i] = {"seq": _refseq, "taxa": _taxa}
                taxa_DB_filesize2 = sys.getsizeof(taxa_DB)

                # Check if the updated taxa DB size is smaller than the previous one and ask for confirmation to update.
                print(taxa_DB_filesize2)
                print(taxa_DB_filesize1)

                if taxa_DB_filesize1 > taxa_DB_filesize2:
                    print("Caution: Updated taxa_DB is smaller size of previous one")
                    choice = input("Are you sure you want to update this? [yes/NO]: ").lower()
                    if choice in ['y', 'ye', 'yes']:
                        print("Updating")
                    elif choice in ['n', 'no']:
                        print("Suspend updating")
                        taxa_DB_path_fail = taxa_DB_path.split(".pkl")
                        with open(taxa_DB_path_fail, "wb") as f:
                            pickle.dump(taxa_DB, f)
                        return

                # Save the updated taxa database to the file.
                with open(taxa_DB_name, "wb") as f:
                    pickle.dump(taxa_DB, f)
            else:
                print("taxa_DB was not updated.")

        def rename_table_DBid():
            # Load the taxa database and the renamed composition table.
            with open(taxa_DB_name, "rb") as f:
                taxa_DB = pickle.load(f)
            comp_df = pd.read_csv(f"{out_dir}{prefix}_rename.csv", index_col = 0)

            # Filter the composition table based on DB IDs present in the taxa database.
            bool_list = [i for i in comp_df.columns if i in taxa_DB.keys()]
            comp_df = comp_df[bool_list]
            bool_list = [i for i in comp_df.columns if taxa_DB[i]["taxa"]["genus"] != "Gifu"]
            comp_df = comp_df[bool_list]
            bool_list = [i for i in comp_df.columns if taxa_DB[i]["taxa"]["family"] != "Chloroplast"]
            comp_df = comp_df[bool_list]

            # Save the filtered composition table as a CSV file.
            comp_df.to_csv(f"{out_dir}{prefix}_microbiome.csv")

            def taxa_aggregate(taxa):
                # Group the composition table by the specified taxa level and sum the counts.
                comp_df_rev = comp_df.T
                taxa_list = list()
                for i in comp_df_rev.index:
                    taxa_list.append(taxa_DB[i]["taxa"][taxa])
                comp_df_rev[taxa] = taxa_list
                comp_df_rev = comp_df_rev.groupby(taxa).sum()
                comp_df_rev = comp_df_rev.T
                comp_df_rev.to_csv(f"{out_dir}{prefix}_microbiome_{taxa}.csv")

            taxa_ref = ['Strain_information', 'domain', 'phylum', 'class', 'order', 'family', 'genus']
            for i in taxa_ref:
                taxa_aggregate(i)


        def bdid2taxa():
            # Load the taxa database and create a DataFrame from it.
            with open(taxa_DB_name, "rb") as f:
                taxa_DB = pickle.load(f)
            temp_dict = dict()
            for i in taxa_DB.keys():
                temp_dict[i] = taxa_DB[i]["taxa"]
            taxa_df = pd.DataFrame.from_dict(temp_dict, orient = "index")
            taxa_df.to_csv(f"{out_dir}{prefix}_taxa.csv")


        # Call the sequence ID update function.
        seqIDupdate()

        # Check if BLAST is required based on the configuration.
        if need_blast == "True":
            blast()
            maketaxadb()
            rename_table_DBid()
            bdid2taxa()
        print("MAUI pipeline was finished.")































#
