""" For a given GOI, returns the specific AA level mutations found
	On a per-cell basis. """

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


def get_filenames_test():
	""" get file names given path """
	files = []
	for file in os.listdir(cwd + "artificalVCF/"):
		if file.endswith(".vcf"):
			fullPath = cwd + 'artificalVCF/' + file 
			files.append(fullPath)
    
	return files



def get_filenames():
	""" get file names given path """
	files = []
	for file in os.listdir(cwd + "scVCF_filtered_subset/"):
		if file.endswith(".vcf"):
			fullPath = cwd + 'scVCF_filtered_subset/' + file 
			files.append(fullPath)
    
	return files



def get_overlap(a, b):
	""" return the len of overlap between two regions """
	q_start = int(a[0])
	q_end = int(a[1])
	r_start = int(b[0])
	r_end = int(b[1])

	if q_start == q_end:
		q_end += 1 		# otherwise it doesnt work for SNPs

	ret = max(0, min(q_end, r_end) - max(q_start, r_start))
	return(ret)



def get_laud_db():
    """ returns the COSMIC database after lung and fathmm filter """
    pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
    database_filter = database.iloc[pSiteList]
    keepRows = database_filter['FATHMM score'] >= 0.7
    db_fathmm_filter = database_filter[keepRows]
    db_fathmm_filter = db_fathmm_filter.reset_index(drop=True)

    return db_fathmm_filter



def write_csv(dictObj, outFile):
	""" writes dict to csv"""
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])



def get_genome_pos_str(sample):
	""" returns a genome position string to match against the ones found in COSMIC """

	chr = sample[0]
	chr = chr.replace("chr", "")
	pos = int(sample[1])
	ref = str(sample[3])
	alt = str(sample[4])

	genomePos = chr + ":" + str(pos) + '-'
	altSplt = alt.split(',')
	add = len(max(altSplt , key = len))

	if add == 1: # for SNPs, start actually == end
		add = 0

	end = pos + add
	genomePos = genomePos + str(end)

	ret = [genomePos, ref, alt]
	return(ret)



def breakdown_by_mutation_type(currSet):
	""" breaks down the mutations set into two categories: SNP and indel
		returns seperate lists for each """

	SNP_list = []
	indel_list = []

	for elm in currSet:
		if len(elm[1]) > 1 or len(elm[2]) > 1: # indel
			indel_list.append(elm)
		else: # SNP
			SNP_list.append(elm)

	retSet = []
	retSet.append(SNP_list)
	retSet.append(indel_list)

	return retSet



def get_mutation_aa(d, chr):
	""" given a dict of {cell, list(genomePos)}, returns a dict of 
		{cell, list(mutation.AA)} 

		THIS NEEDS WORK!! """

	print('AA searching...')
	newDict = {}

	for k in d:
		valuesList = d.get(k) # can now handle values with multiple entries
		newValues = []

		for entry in valuesList:
			testSplit = entry.split('-') # if its a SNP it wont have '-' at all	

			### CASE 1 -- SNP
			if len(testSplit) == 1:
				chrStr = chr + ':' + entry + '-' + entry
				filter_df = database_laud[database_laud["Mutation genome position"].str.contains(chrStr)==True]
				currMuts = filter_df['Mutation AA']

				for item in currMuts:		# really shouldnt have a for loop here
					item = item.replace("p.", "")

				newValues.append(item) 		# effectively just taking the last item in the list
			
			### CASE 2 -- INDEL 
			else:
				chrStr = chr + ':' + entry
				filter_df = database_laud[database_laud["Mutation genome position"].str.contains(chrStr)==True]
				currMuts = filter_df['Mutation AA']

				for item in currMuts:		# really shouldnt have a for loop here
					item = item.replace("p.", "")

				newValues.append(item)		# effectively just taking the last item in the list

		newDict.update({k : newValues})

	return newDict



def get_shared_ROIs(currSet, genomePos_laud_db_, SNP_bool):
	""" given a set of genome position strings, searches for the ones
		that are in the COSMIC database """
	if SNP_bool: # exact string match
		curr_df = pd.DataFrame(currSet, columns=['posStr', 'ref', 'alt'])
		overlap = set(curr_df['posStr']).intersection(set(genomePos_laud_db_))
		keep = curr_df.where(curr_df['posStr'].isin(overlap))
		keep = keep.dropna()

		ret = keep.values.tolist() # convert entire df to list

	else: # set intersect -- NOT SURE HOW TO DO THIS YET
		ret = [['0:0-0', 'N', 'N']] # dummy

	return(ret)



def hit_search_genome_coords(sample, *args):
	""" given a list of shared entries between an individual cell's VCF 
		and the COSMIC LAUD db, searches for hits to the GOI """
	cell_ = args[0]
	queryChrom_ = args[1]
	lPosQuery_ = args[2]
	rPosQuery_ = args[3]
	match = ""

	posStr = sample[0]

	chrom = posStr.split(':')[0]
	start = posStr.split(':')[1].split('-')[0]
	end = posStr.split(':')[1].split('-')[1]
	ref = posStr[1]
	alt = posStr[2]

	if chrom == queryChrom_:
		overlap = get_overlap([start, end], [lPosQuery_, rPosQuery_])
	
		if overlap != 0: # curr sample is in the gene of interest!
			print(overlap)
			match = posStr
	
	return match



def driver(fileNames, chrom, pos1, pos2):
	""" creates dict with genome coords for hits to specific GOI """
	print('getting coords to GOI hits')

	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])
	cells_dict_GOI_coords = {}
	queryChrom = chrom
	lPosQuery = pos1
	rPosQuery = pos2

	for f in fileNames:
		cell = f.replace("/home/ubuntu/cerebra/cerebra/wrkdir/scVCF_filtered_subset/", "")
		cell = cell.replace(".vcf", "")	
	
		df = VCF.dataframe(f)
		genomePos_query = df.apply(get_genome_pos_str, axis=1) # apply function for every row in df
		genomePos_query_breakdown = breakdown_by_mutation_type(genomePos_query)
		
		SNPs = genomePos_query_breakdown[0]
		indels = genomePos_query_breakdown[1]

		# get the entries shared between curr cells VCF and the LAUD filter set
		#	remember, these are general, and NOT gene specific
		shared_SNPs = get_shared_ROIs(SNPs, genomePos_laud_db, 1)
		shared_indels = get_shared_ROIs(indels, genomePos_laud_db, 0)
		
		all_shared_regions = shared_SNPs + shared_indels

		shared_pd = pd.Series(all_shared_regions)
		matches = shared_pd.apply(hit_search_genome_coords, args=(cell, queryChrom, lPosQuery, rPosQuery))

		# delete empty dict keys
		for k in matches.keys():
			try:
				if len(matches[k])<1:
					del matches[k]
			except: pass

		cells_dict_GOI_coords.update({cell : list(matches.values)})

	return cells_dict_GOI_coords



""" get cmdline input """
@click.command()
@click.option('--test', default = False)
@click.option('--chrom', default = 7, prompt='chromosome', required=True, type=str)
@click.option('--start', default = 55152337, prompt='start position', required=True, type=int)
@click.option('--end', default = 55207337, prompt='end position', required=True, type=int)
@click.option('--outprefix', default = 'sampleOut', prompt='prefix to use for outfile', required=True, type=str)
@click.option('--wrkdir', default = '/home/ubuntu/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)
 


def get_specific_mutations(test, chrom, start, end, outprefix, wrkdir):
	""" for a specific gene of interest, get the complete set of amino acid level mutations
		for each cell in dataset """
	global database
	global database_laud
	global cwd
	global test_bool

	cwd = wrkdir
	test_bool = test
	chrom = str(chrom)
	start = str(start)
	end = str(end)

	print('setting up COSMIC database...')
	database = pd.read_csv(cwd + "CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = get_laud_db()

	if test_bool:
		fNames = get_filenames_test()
	else:
		fNames = get_filenames()

	goiDict = driver(fNames, chrom, start, end) # get genome coords
	#print(goiDict)
	print("GOI search done!")

	goiDict_AA = get_mutation_aa(goiDict, chrom)
	print('AA search done')	

	if test_bool:
		write_csv(goiDict, cwd + 'test/specific_mutations/' + outprefix + '.csv')
		write_csv(goiDict_AA, cwd + 'test/specific_mutations/' + outprefix + '_AA.csv')
	else:
		write_csv(goiDict, cwd + outprefix + '.csv')
		write_csv(goiDict_AA, cwd + outprefix + '_AA.csv')
