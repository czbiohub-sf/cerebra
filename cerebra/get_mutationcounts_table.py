""" 
This program generates a gene/cell counts table for the mutations
found in a given set of vcf files
"""

import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import csv
import pandas as pd
import sys
import multiprocessing as mp
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_germline_filtered_cells_list():
	""" returns a list of cells that have already been germline filterd """
	filterDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/filteredOut/'
	filterDir_list = os.listdir(filterDir)

	filteredCells = []
	for f in filterDir_list:
		cell = f.strip('_unique.vcf')
		filteredCells.append(cell)

	return filteredCells



def get_file_names():
	""" get file names based on specified path """
	files = []
	for file in os.listdir("vcf_test/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("vcf_test/", file))
			files.append(fullPath)
    
	return files



def get_laud_db():
	""" returns the COSMIC database after the lung adeno filter """
	print('setting up LAUD filtered database...')
	pHistList = database.index[database['Primary histology'] == 'carcinoma'].tolist()
	pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
	shared = list(set(pHistList) & set(pSiteList))
	database_filter = database.iloc[shared]
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
		genomePos = 'ERROR'

	return(genomePos)



def get_gene_name(posString):
	""" return the gene name from a given genome position string
	   (ie. '1:21890111-21890111'), by querying the hg38-plus.gtf """

	chrom = posString.split(':')[0] # work on posString
	posString_remove = posString.split(':')[1]
	lPosition = posString_remove.split('-')[0] 
	rPosition = posString_remove.split('-')[1] 

	# work on hg38_gtf
	chromStr = 'chr' + str(chrom)
	hg38_gtf_filt = hg38_gtf.where(hg38_gtf[0] == chromStr).dropna()
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[3] <= int(lPosition)).dropna() # lPos good
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[4] >= int(rPosition)).dropna() # rPos good
	
	try:
		returnStr = str(hg38_gtf_filt.iloc[0][8])	# keep just the gene name / meta data col
		returnStr = returnStr.split(';')[1]
		returnStr = returnStr.strip(' gene_name')
		returnStr = returnStr.strip(' ')
		returnStr = returnStr.strip('"')
	except IndexError:
		returnStr = ''
	
	return returnStr



def get_genecell_mut_counts(f):
	""" creates a dictionary obj where each key is a cell and each value
	is a list of the genes we found mutations for in that cell """
	tup = [] 

	cell = f.replace("vcf_test/", "")
	cell = cell.replace(".vcf", "")
	print(cell) # to see where we are
	
	df = VCF.dataframe(f)
	genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df

	items = set(genomePos_query) # genomePos_query (potentially) has dups

	# COSMIC filter
	shared = [i for i in genomePos_laud_db if i in items] # retains dups

	shared_series = pd.Series(shared)
	sharedGeneNames = shared_series.apply(getGeneName)
	tup = [cell, sharedGeneNames]

	return(tup)



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

	return(df)



""" get cmdline input """
@click.command()
@click.option('--COSMIC_db', default = '', prompt='s3 path to COSMIC database', required=True, type=str)
@click.option('--vcf_path', default = '', prompt='s3 path to scVCFs', required=True, type=str)
@click.option('--outpath', default = 's3://lincoln.harris-work/filteredOut/', prompt='s3 path to where output table should be pushed', required=True, type=str)



def get_mutationcounts_table(database, vcf_path, outpath):
	""" driver function """
	global database
	global database_laud
	global hg38_gtf
	global genomePos_laud_db
	global germlineFilterCells

	database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = getLAUD_db()
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])
	hg38_gtf = pd.read_csv('hg38-plus.gtf', delimiter = '\t', header = None)
	fNames = getFileNames()

	germlineFilterCells = getGermlineFilteredCellsList()

	print('creating pool')

	p = mp.Pool(processes=14)

	try:
		cells_list = p.map(getGeneCellMutCounts, fNames, chunksize=1) # default chunksize=1
	finally:
		p.close()
		p.join()

	# convert to dictionary
	cells_dict = {}

	for item in cells_list:
    	cells_dict.update({item[0]:item[1]})

	print('writing file')

	filterDict_pd = pd.DataFrame.from_dict(cells_dict, orient="index") # orient refers to row/col orientation 
	filterDict_format = formatDataFrame(filterDict_pd)
	filterDict_format.to_csv("geneCellMutationCounts.csv")
