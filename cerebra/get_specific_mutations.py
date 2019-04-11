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


def get_filenames():
	""" get file names given path """
	files = []
	for file in os.listdir("wrkdir/scVCF_filtered_all/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("wrkdir/scVCF_filtered_all/", file))
			files.append(fullPath)
    
	return files



def get_genome_pos(sample):
	""" returns a genome position string to match against the ones found in COSMIC """
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



def get_laud_db():
	""" returns the COSMIC database after lung adeno filter """
	print('setting up LAUD filtered database...')
	pHistList = database.index[database['Primary histology'] == 'carcinoma'].tolist()
	pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
	shared = list(set(pHistList) & set(pSiteList))
	database_filter = database.iloc[shared]
	return database_filter



def hit_search(sample):
	""" performs the actual search """
	match = 0
	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on `-`
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1
			# rPosQuery and lPosQuery are GLOBALs
			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left position good
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right position good
					match = 1
		except IndexError:
			print('index error')

	return match



def hit_search_coords(sample, *args):
	""" given a list of shared entries between an individual cell's VCF 
		and the COSMIC LAUD db, searches for hits to the GOI """
	cell_ = args[0]
	match = ""

	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on `-`
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1
			# keep in mind rPosQuery and lPosQuery are GLOBALs
			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left pos GOI match
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right pos GOI match
					if lPosCurr == rPosCurr: # SNP
						match = lPosCurr
					else: 		# found an indel!!
						match = lPosCurr + '-' + rPosCurr
						#print(cell_)

		except IndexError:
			print('index error')

	return match



def get_goi_hits(fileNames, chrom, pos1, pos2):
	""" creates dict with hits to a specific GOI """
	print('getting hits to GOI')

	global queryChrom, lPosQuery, rPosQuery # dont like this
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	cells_dict_GOI = {}
	queryChrom = chrom
	lPosQuery = pos1
	rPosQuery = pos2

	for f in fileNames:
		numMatches = 0
		cell = f.replace("vcf_germline_filter/", "")
		cell = cell.replace(".vcf", "")	

		df = VCF.dataframe(f)
		genomePos_query = df.apply(get_genome_pos, axis=1) # apply function for every row in df
		
		shared = list(set(genomePos_query) & set(genomePos_laud_db)) # get the LAUD filter set
		shared1 = pd.Series(shared) # what if i convert this guy to a pandas object? 

		numMatches = shared1.apply(hit_search_func) # another apply call 

		cells_dict_GOI.update({cell : sum(numMatches)})
	
	return cells_dict_GOI



def get_goi_hits_coords(fileNames, chrom, pos1, pos2):
	""" creates dict with genome coords for hits to specific GOI """
	print('getting coords to GOI hits')

	global queryChrom, lPosQuery, rPosQuery # dont like this
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])
	cells_dict_GOI_coords = {}
	queryChrom = chrom
	lPosQuery = pos1
	rPosQuery = pos2

	for f in fileNames:
		numMatches = 0
		cell = f.replace("vcf_germline_filter/", "")
		cell = cell.replace(".vcf", "")	

		df = VCF.dataframe(f)
		genomePos_query = df.apply(get_genome_pos, axis=1) # apply function for every row in df
		# get the entries shared between curr cells VCF and the LAUD filter set
		#	remember, these are general, and NOT gene specific
		genomePos_query_expand = expandSet(set(genomePos_query))

		shared = list(set(genomePos_query_expand) & set(genomePos_laud_db))
		shared1 = pd.Series(shared) # convert to pandas obj
		matches = shared1.apply(hit_search_coords, args=(cell,)) # another apply call 

		# delete empty dict keys
		for k in matches.keys():
			try:
				if len(matches[k])<1:
					del matches[k]
			except: pass

		cells_dict_GOI_coords.update({cell : list(matches.values)})

	return cells_dict_GOI_coords



def get_mutation_aa(d, chr):
	""" given a dict of {cell, list(genomePos)}, returns a dict of 
		{cell, list(mutation.AA)} """
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



def expand_set(mySet):
	""" pass in a set of genome coords, and it will 'expand' the indels
		within the set by adding +/- 3 bp copies for each one """
	returnSet = []

	for entry in mySet:
		l0 = []
		l1 = []
		try:
			sub0 = entry.split('-')[0] # split on `-`
			sub1 = entry.split('-')[1] # this guy is good
			sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom
			chrom = sub0.split(':')[0]
		
			if sub00 != sub1: # got an indel 
				sub00_1 = int(sub00) + 1
				sub00_2 = int(sub00) + 2
				sub00_3 = int(sub00) + 3
				sub00_4 = int(sub00) - 1
				sub00_5 = int(sub00) - 2
				sub00_6 = int(sub00) - 3

				l0.extend((sub00_1, sub00_2, sub00_3, sub00_4, sub00_5, sub00_6))
				
				try:
					sub1_1 = int(sub1) + 1
					sub1_2 = int(sub1) + 2
					sub1_3 = int(sub1) + 3
					sub1_4 = int(sub1) - 1
					sub1_5 = int(sub1) - 2
					sub1_6 = int(sub1) - 3

					l1.extend((sub1_1, sub1_2, sub1_3, sub1_4, sub1_5, sub1_6))
				
				except ValueError:
					continue

				coord_combos = list(itertools.product(l0, l1))
				for pair in coord_combos:
					toAdd = chrom + ':' + str(pair[0]) + '-' + str(pair[1])
					returnSet.append(toAdd)

			else:
				returnSet.append(entry)
		
		except IndexError:
			continue

	return returnSet



def writeCSV(dictObj, outFile):
	""" writes dict to csv"""
	print('writing csv')
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])



""" get cmdline input """
@click.command()
@click.option('--chrom', default = 7, prompt='chromosome', required=True, type=int)
@click.option('--start', default = 55152337, prompt='start position', required=True, type=int)
@click.option('--end', default = 55207337, prompt='end position', required=True, type=int)
@click.option('--outfile', default = 'sampleOut.csv', prompt='name of output file', required=True, type=str)
 


def get_specific_mutations(chrom, start, end, outfile):
	""" driver func """
	global database
	global database_laud

	print('setting up COSMIC database...')
	database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = get_laud_db()
	#database_laud.to_csv('database_laud.csv') # test write
	fNames = get_filenames()

	chromo = sys.argv[2]
	position1 = sys.argv[3]
	position2 = sys.argv[4]	

	#goiDict = getGOIHits(fNames, chromo, position1, position2) # standard call - get raw counts
	goiDict = get_goi_hits_coords(fNames, chromo, position1, position2) # get genome coords
	print("GOI search done!")

	outFilePref = sys.argv[5]
	write_csv(goiDict, outFilePref + '.csv')

	goiDict_AA = get_mutation_aa(goiDict, chromo)
	print('AA search done')
	write_csv(goiDict_AA, outFilePref + '_AA.csv')

