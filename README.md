# Genome_Quebec_Interactive_Sanger_Analysis

#Working directory setup requirements:

#Minimum requirements: a folder named Ab1 with .ab1 type sanger sequencing files in it
#Optional folders:
#plate_info with one or more genome quebec plate data .csv files
#this file will result in ab1 files being automatically filtered based on Genome quebec quality control results

#References with one or more .txt files of reference genes in fasta format
#This folder will enable optional alignment of reference sequences with the imported sanger seq results

#metadata with one or more .csv files of sample metadata
#This folder will enable optional output of combined metadat/sequencing result csvs
#NOTE: the metadata .csvs must include a column with your full sample names (the same names that were sent in the 
#genome quebec request excel sheet)
