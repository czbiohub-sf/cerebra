""" 
This program generates a gene/cell counts table for the mutations
found in a given set of vcf files
"""

import numpy as np
import vcf
import os
import csv
import pandas as pd
import sys
import multiprocessing as mp
import warnings
import click
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_filenames_test():
	""" get file names, for the test condition"""
	files = []
	for file in os.listdir(cwd + "artificalVCF/"):
		PATH = cwd + 'artificalVCF/' + file
		files.append(PATH)

	return files



def get_filenames():
	""" get file names based on specified path """
	files = []
	for file in os.listdir(cwd + "scVCF_filtered_subset/"):
		PATH = cwd + 'scVCF_filtered_subset/' + file
		files.append(PATH)

	return files



def get_laud_db():
	""" returns the COSMIC database after lung adeno filter """
	print('setting up LAUD filtered database...')
	#pHistList = database.index[database['Primary histology'] == 'carcinoma'].tolist()
	pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
	#shared = list(set(pHistList) & set(pSiteList))
	database_filter = database.iloc[pSiteList]

	return database_filter



def get_genome_pos(sample):
	""" returns the genome position string that will match against the 
		ones in COSMIC db  """
	try:
		chr = str(sample[0])
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
		else: # multibase-for-multibase substitution
			secondPos = '1'
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	except:
		# NOTE: Should probably throw here.
		genomePos = 'ERROR'

	return genomePos

def parse_genome_pos(pos_str):
	chrom, pos_range = pos_str.split(':')[0:2]
	start_pos, end_pos = pos_range.split('-')[0:2]

	return (chrom, start_pos, end_pos)

def make_genome_pos(record):
	# Although including `chr` in the CHR column constitutes malformed VCF, it
	# may be present, so it should be removed.
	CHROM = record.CHROM.replace("chr", "")
	POS = record.POS
	ref_len = len(record.REF)
	alt_len = len(record.ALT)

	if ref_len == 1 and alt_len == 1:
		return (CHROM, POS, POS)
	elif ref_len > 1 and alt_len == 1:
		return (CHROM, POS, POS + ref_len)
	elif alt_len > 1 and ref_len == 1:
		return (CHROM, POS, POS + alt_len)
	else: # multibase-for-multibase substitution
		return (CHROM, POS, 1)

def find_gene_name(genome_pos):
	""" return the gene name from a given genome position string
	   (ie. '1:21890111-21890111'), by querying the hg38-plus.gtf """

	chrom, pos_start, pos_end = genome_pos

	# work on hg38_gtf
	chromStr = 'chr' + str(chrom)
	hg38_gtf_filt = hg38_gtf.where(hg38_gtf[0] == chromStr).dropna()
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[3] <= int(pos_start)).dropna() # lPos good
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[4] >= int(pos_end)).dropna() # rPos good
	
	try:
		returnStr = str(hg38_gtf_filt.iloc[0][8])	# keep just the gene name / metadata col
		returnStr = returnStr.split(';')[1]
		returnStr = returnStr.strip(' gene_name')
		returnStr = returnStr.strip(' ')
		returnStr = returnStr.strip('"')
	except IndexError:
		returnStr = ''

	return returnStr



def find_cell_mut_gene_names(filename):
	""" creates a dictionary obj where each key is a cell and each value
		is a list of the genes we found mutations for in that cell """
	tup = []


	cell = filename.replace(cwd, "")
	cell = cell.replace('scVCF_filtered_all/', "")
	cell = cell.replace(".vcf", "")

	# PyVCF documentation claims that it automatically infers compression type
	# from the file extension.
	vcf_reader = vcf.Reader(filename=filename)

	filtered_gene_names = []

	for record in vcf_reader:
		genome_pos = make_genome_pos(record)
		# TODO: Filter out duplicates?
		# And is it better to filter out positional duplicates or gene name
		# duplicates?

		# No COSMIC filter in test mode
		if not test_bool:
			# Skip this gene if it isn't in the laud_db
			if genome_pos in genome_pos_laud_db:
				continue
		
		gene_name = find_gene_name(genome_pos)

		filtered_gene_names.append(gene_name)

	filtered_series = pd.Series(filtered_gene_names)

	tup = (cell, filtered_series)

	return tup



def format_dataframe(raw_df):
	""" creates the cell/mutation counts table from the raw output that 
		get_gene_cell_muts_counts provides """
	cellNames = list(raw_df.index)

	genesList = []
	for i in range(0, raw_df.shape[0]):
		currList = list(raw_df.iloc[i].unique()) # unique genes for curr_cell 

		for elm in currList:
			if elm not in genesList:
				genesList.append(elm)
	
	genesList1 = pd.Series(genesList)

	df = pd.DataFrame(columns=genesList1, index=cellNames) # initialize blank dataframe
	for col in df.columns: # set everybody to zero
		df[col] = 0

	for i in range(0,raw_df.shape[0]):
		currCell = raw_df.index[i]
		currRow = raw_df.iloc[i]

		for currGene in currRow:
			df[currGene][currCell] += 1

	return df



""" get cmdline input """
@click.command()
@click.option('--nthread', default = 64, prompt='number of threads', required=True, type=int)
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/home/ubuntu/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def get_mutationcounts_table(nthread, test, wrkdir):
	""" generate a cell x gene mutation counts table from a set of germline filtered vcfs """
	global database
	global database_laud
	global hg38_gtf
	global genome_pos_laud_db
	global cwd
	global test_bool

	cwd = wrkdir
	test_bool = test

	database = pd.read_csv(cwd + "CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = get_laud_db()
	genome_pos_laud_db = set(map(parse_genome_pos, database_laud['Mutation genome position']))
	hg38_gtf = pd.read_csv(cwd + 'hg38-plus.gtf', delimiter = '\t', header = None)

	if test:
		fNames = get_filenames_test()
	else:
		fNames = get_filenames()
	
	print('creating pool')
	p = mp.Pool(processes=nthread)
	print('running...')

	try:
		cell_genes_pairs = p.map(find_cell_mut_gene_names, fNames, chunksize=1) # default chunksize=1
	finally:
		p.close()
		p.join()

	# convert to dictionary
	cells_dict = {}
	naSeries = pd.Series([np.nan])

	for cell, gene_names in cell_genes_pairs:
		if len(gene_names.index) == 0:
			toAdd = {cell:naSeries}
		else:
			toAdd = {cell:gene_names}
		cells_dict.update(toAdd)

	print('writing file')

	filterDict_pd = pd.DataFrame.from_dict(cells_dict, orient="index") # orient refers to row/col orientation 
	filterDict_format = format_dataframe(filterDict_pd)

	filterDict_format.to_csv(cwd + "intermediate.csv")
	intermediate = pd.read_csv(cwd + 'intermediate.csv')

	# rename 0th col
	intermediate.rename(columns={'Unnamed: 0' :'cellName'}, inplace=True)

	cols = intermediate.columns
	dropList = []

	for col in cols:
		if 'Unnamed' in col:
			dropList.append(col)

	# want to drop the cols that contain 'Unnamed'
	intermediate = intermediate.drop(dropList, axis=1)

	if test_bool:
		intermediate.to_csv(cwd + "test/mutationcounts_table/geneCellMutationCounts_artifical.csv", index=False)	
	else:
		intermediate.to_csv(cwd + "geneCellMutationCounts_noLaudFilter.csv", index=False)

	cmd = 'rm ' + cwd + 'intermediate.csv'
	os.system(cmd)
