#!/usr/bin/env python
import os
import pandas as pd
from glob import glob
import re
import subprocess

# Identify folders from list of folders
def matchFolders(x, dir_list = None):
    for folder in dir_list:
        # Build regex expression
        regex_exprs1 = re.compile("/" + x + "/$")
        regex_exprs2 = re.compile("/" + x + "_\w*/$")
        if re.search(regex_exprs1, folder):
            return(folder)
        elif re.search(regex_exprs2, folder):
            return(folder)
        else:
            return(None)

def get_scrnaseq_dir_list(scrnaseq_dir, dir_list = None, pools = None):
    try:
        scrnaseq_glob_exprs = os.path.join(scrnaseq_dir,  "*/")
        query_directories = glob(scrnaseq_glob_exprs)
        scrnaseq_filelist = [matchFolders(pool, dir_list = query_directories) for pool in pools]
        scrnaseq_filedict = dict(zip(pools, scrnaseq_filelist))
        scrnaseq_libs = pd.Series(scrnaseq_filedict, name="scRNAseq_Directories")
        return(scrnaseq_filelist, scrnaseq_libs)
    except TypeError as error:
        print(error)
        raise SystemExit("Could not find a scRNA-seq directory for all of the pools in your pool list. Please check that they are spelled correctly and you do not have any additional pool names that are not in {}  ".format(individual_list_dir))

def get_barcodes_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("barcodes.tsv", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))

### Get the barcode files for each pool     
def get_barcodes_dir(scrnaseq_filelist, pools = None):    
    try:
        barcode_filelist = [get_barcodes_files(pool) for pool in scrnaseq_filelist]
        barcode_filedict = dict(zip(pools, barcode_filelist))
        barcode_libs = pd.Series(barcode_filedict, name="Barcode_Files")
        return(barcode_libs)
    except Exception as error:
        print(error)
        raise SystemExit("Could not find a barcode file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'barcodes.tsv' within the name.")
    
def get_bam_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if f.endswith(".bam")]:
            return(os.path.join(dirpath, filename))

def get_bam_dirs(scrnaseq_filelist, pools = None):
    try:
        bam_filelist = [get_bam_files(pool) for pool in scrnaseq_filelist]
        bam_filedict = dict(zip(pools, bam_filelist))
        bamlibs = pd.Series(bam_filedict, name="Bam_Files")
        return(bamlibs)
    except Exception as error:
        print(error)
        raise SystemExit("Could not find a bam file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain '.bam' within the name.")
    
def get_matrix_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("matrix.mtx", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))

def get_matrix_dirs(scrnaseq_filelist, pools = None):
    try:
        matrix_filelist = [get_matrix_files(pool) for pool in scrnaseq_filelist]
        matrix_filedict = dict(zip(pools, matrix_filelist))
        matrix_libs = pd.Series(matrix_filedict, name="Matrix_Files")
        matrix_dirlist = [os.path.dirname(filename) for filename in matrix_filelist]
        matrix_dirdict = dict(zip(pools, matrix_dirlist))
        matrix_dir_libs = pd.Series(matrix_dirdict, name="Matrix_Directories")
        return(matrix_dir_libs, matrix_libs)
    except Exception as error:
        print(error)
        raise SystemExit("Could not find a matrix file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'matrix.mtx' within the name.")
        

def get_individual_files(x, individual_dir = None):
    for filename in individual_dir:
        if re.search(r'^' + x + "$", filename):
            return(filename)
        elif re.search(r'^' + x + '\D', filename):
            return(filename)
        elif re.search(x + "$", filename):
            return(filename)
        elif re.search(x + '\D', filename):
            return(filename)

def get_individual_dirs(individual_list_dir, pools = None):
    individual_dirlist = os.listdir(individual_list_dir)
    try:
        individual_filelist = [os.path.join(individual_list_dir, get_individual_files(pool, individual_dir = individual_dirlist)) for pool in pools]
        individual_filedict = dict(zip(pools, individual_filelist))
        individual_libs = pd.Series(individual_filedict, name="Individuals_Files")
        return(individual_libs)
    except Exception as error:
        print(error)
        raise SystemExit("Could not find a files of individuals in {}. Please check that they exist somewhere in this directory and contain the pool names within the name of the file.".format(individual_list_dir))
    
def get_scrnaseq_dirs(config):
    # Extract variables from configuration file for use within the rest of the pipeline
    input_dict = config["inputs"]
    output_dict = config["outputs"]
    ref_dict = config["refs"]

    # General variables used by the rest of the pipeline
    scrnaseq_dir = input_dict["scRNAseq_dir"]
    individual_list_dir = input_dict["individual_list_dir"]

    ### Check that all the directories exist and can be accessed
    if not os.path.exists(scrnaseq_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(scrnaseq_dir))
    if not os.path.exists(individual_list_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(individual_list_dir))

    # Read in samplesheet from the configfile specified from the command line
    samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")

    # Expect first colunn to be pools
    pools = samples.iloc[:, 0]

    # Match pools to scrna seq directories to make a list of each scRNA-seq dir
    scrna_seq_dirlist = [d for d in os.listdir(scrnaseq_dir) if os.path.isdir(os.path.join(scrnaseq_dir, d))]
    scrnaseq_filelist, scrnaseq_libs = get_scrnaseq_dir_list(scrnaseq_dir, dir_list=scrna_seq_dirlist, pools = pools)
    
    # Get files from scrnaseq directories
    barcode_libs = get_barcodes_dir(scrnaseq_filelist, pools = pools)
    
    # Get the bam files for each pool
    bamlibs = get_bam_dirs(scrnaseq_filelist, pools = pools)

    ### Get the matrix files for each pool
    matrix_dir_libs, matrix_libs = get_matrix_dirs(scrnaseq_filelist, pools = pools)

    ### Get the matrix individual list file for each pool
    individual_libs = get_individual_dirs(individual_list_dir, pools = pools)
    dataframe = pd.concat([scrnaseq_libs, barcode_libs, bamlibs, matrix_libs, matrix_dir_libs, individual_libs], axis=1)

    # Check if any files are missing
    if dataframe.isnull().values.any():
        missing_rows = (dataframe.isnull()).any(axis = 1)
        missing_samples = missing_rows.index.to_list()
        missing_sample_names = (",").join(missing_samples)
        raise SystemExit("The following samples are missing data: {missing_sample_names}".format(missing_sample_names = missing_sample_names))
    else:
        return(dataframe)

