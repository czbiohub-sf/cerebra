""" Gets the total number of mutations, on a per-cell basis.  
	Two run modes -- raw and filter (COSMIC, LAUD annotation) """

import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import csv
import pandas as pd
import sys
import itertools
import warnings
import click 
warnings.simplefilter(action='ignore', category=FutureWarning)



def get_filenames():
	""" gets file names given path """
	files = []
	for file in os.listdir("wrkdir/scVCF_filtered_all/"):
		fullPath = (os.path.join("wrkdir/scVCF_filtered_all/", file))
		files.append(fullPath)
    
	return files



def get_raw_counts(fileNames):
	""" returns dict obj with raw counts for GATK hits w/in set of vcfs"""
	print('getting raw counts...')
	cells_dict = {}

	for f in fileNames:
		cell = f.replace("wrkdir/scVCF_filtered_all/", "")
		cell = cell.replace(".vcf", "")
    
		df = VCF.dataframe(f)
		unique = len(np.unique(df.POS))
    
		cells_dict.update({cell : unique})
	print('finished!')
	return cells_dict



def get_genome_pos(sample):
	""" returns genome position string that will match against the one in COSMIC"""
	try:
		chr = sample[0]
		chr = chr.replace("chr", "")
		pos = int(sample[1])
		ref = str(sample[3])
		alt = str(sample[4])
	
		if (len(ref) == 1) & (len(alt) == 1): # most basic case
			secondPos = pos
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		elif (len(ref) > 1) & (len(alt) == 1):
			secondPos = pos + len(ref)
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		elif (len(alt) > 1) & (len(ref) == 1):
			secondPos = pos + len(alt)
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		else: # BOTH > 1 .... not sure what to do here. does this actually happen? 
			secondPos = 'dummy'
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	except:
		genomePos = 'chr0:0-0'

	return(genomePos)



def get_filter_counts_basic(fileNames):
	""" returns dict with COSMIC filtered GATK hits """
	print('getting filter counts basic...')
	cells_dict_filter = {}
	genomePos_db = pd.Series(database['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("wrkdir/scVCF_filtered_all/", "")
		cell = cell.replace(".vcf", "")
		df = VCF.dataframe(f)

		genomePos_query = df.apply(get_genome_pos, axis=1)
    
		shared = list(set(genomePos_query) & set(genomePos_db))
		cells_dict_filter.update({cell : len(shared)})
    
	print('finished!')
	return cells_dict_filter



def get_laud_db():
	""" returns the COSMIC database after LAUD filter """
	print('setting up LAUD filtered database...')
	pHistList = database.index[database['Primary histology'] == 'carcinoma'].tolist()
	pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
	shared = list(set(pHistList) & set(pSiteList))
	database_filter = database.iloc[shared]
	return database_filter



def get_filter_counts_laud(fileNames):
	""" returns dict with COSMIC filtered GATK hits """
	print('getting filter counts LAUD...')
	cells_dict_laud = {}
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("wrkdir/scVCF_filtered_all/", "")
		cell = cell.replace(".vcf", "")

		df = VCF.dataframe(f)
		genomePos_query = df.apply(get_genome_pos, axis=1) # apply function for every row in df
    
		shared = list(set(genomePos_query) & set(genomePos_laud_db))
		cells_dict_laud.update({cell : len(shared)})

	print('finished!')
	return cells_dict_laud



def write_csv(dictObj, outFile):
	""" writes dict to csv """
	print('writing csv')
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])



""" get cmdline input """
@click.command()
@click.option('--mode', default = 1, prompt='run mode', required=True, type=int)



def get_mutationalburden(mode):
	""" returns the total number of mutations (mutational burden) for each cell in dataset.
		three run modes: 1) raw counts, 2) COSMIC filter counts, 3) COSMIC lung adenocarcinoma
		filter  """
	global database
	global database_laud

	# raw counts
	if mode == 1:
		fNames = get_filenames()
		rawDict = get_raw_counts(fNames)
		print("raw counts done!")
		print('writing csv')
		write_csv(rawDict, "wrkdir/nonImmune_GATK_hits_raw.csv")

	# filter counts (basic)
	if mode == 2:
		print('setting up COSMIC database...')
		database = pd.read_csv("wrkdir/CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
		fNames = get_filenames()
		filterDict = get_filter_counts_basic(fNames)
		print("filter counts (basic) done!")
		print('writing csv')
		write_csv(filterDict, "wrkdir/nonImmune_GATK_hits_COSMIC_filter.csv")

	# filter counts LAUD
	if mode == 3:
		print('setting up COSMIC database...')
		database = pd.read_csv("wrkdir/CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
		database_laud = get_laud_db()
		fNames = get_filenames()
		filterDict1 = get_filter_counts_laud(fNames) 
		print("filter counts (LAUD) done!")
		print('writing csv')
		write_csv(filterDict1, "wrkdir/nonImmune_GATK_hits_LAUD_filter.csv")
