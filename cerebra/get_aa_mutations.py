""" For a given GOI, returns the specific AA level mutations found
	On a per-cell basis. """

import numpy as np
import os
from . import utils
#import .utils
import csv
import pandas as pd
import sys
import itertools
import warnings
import click 
from tqdm import tqdm
import multiprocessing as mp
warnings.simplefilter(action='ignore', category=FutureWarning)


def count_comments(filename):
	""" Count comment lines (those that start with "#") 
		cribbed from slowko """
	comments = 0
	fn_open = gzip.open if filename.endswith('.gz') else open
	with fn_open(filename) as fh:
		for line in fh:
			if line.startswith('#'):
				comments += 1
			else:
				break
	return comments



def vcf_to_dataframe(filename):
	""" Open a VCF file and return a pandas.DataFrame with
		each INFO field included as a column in the dataframe 
		cribbed from slowko """
	VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20']

	# Count how many comment lines should be skipped.
	comments = count_comments(filename)
	tbl = pd.read_table(filename, compression=None, skiprows=comments,
							names=VCF_HEADER, usecols=range(10))
	
	return(tbl)



def get_filenames():
	""" get file names given path """
	files = []
	for file in os.listdir(cwd + "scVCF_filtered_all/"):
		if file.endswith(".vcf"):
			fullPath = cwd + 'scVCF_filtered_all/' + file 
			files.append(fullPath)
    
	return files



def get_laud_db(gene_):
    """ returns the COSMIC database after lung and fathmm filter  """
    pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
    db_filter = database.iloc[pSiteList]
    #keepRows = db_filter['FATHMM score'] >= 0.7
    #db_fathmm_filter = dbfilter[keepRows]
    #db_fathmm_filter = db_fathmm_filter.reset_index(drop=True)
    keep = db_filter['Gene name'] == gene_
    db_gene = db_filter[keep]
    db_gene = db_gene.reset_index(drop=True)

    if len(db_gene.index) == 0:
    	print('   this gene is not in the cosmic database')
    	print('   maybe it has a different name? ')
    	print('')
    	sys.exit()

    return db_gene
    #return db_fathmm_filter



def write_csv(dictObj, outFile):
	""" writes dict to csv"""
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])



def generate_genome_pos_str(sample):
	""" given vcf record, returns a genome position string"""

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

	if len(ref) > 1: # deletion case
		add = len(ref)

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



def get_corresponding_aa_subs(d):
	""" given a dict of {cell, list(genomePos)}, returns a dict of 
		{cell, list(mutation.AA)} """

	print('finding corresponding amino acid sequences...')
	print(' ')
	newDict = {}

	for cell in d:
		valuesList = d.get(cell) # handles cells with multiple hits
		newValues = []

		for entry in valuesList:
			posStr = entry[0]
			ref = entry[1]
			alt = entry[2]
			nucSub = ref + '>' + alt

			curr_obj = utils.GenomePosition.from_str(posStr)

			if len(ref) > 1 or len(alt) > 1: # indel case -- just return best overlap
				b = cosmic_genome_tree.get_best_overlap(curr_obj)
				if b is not None:
					AA_sub = b["Mutation AA"]
					AA_sub = AA_sub.replace("p.", "")
					newValues.append(AA_sub)
			else: # SNP case -- specific CDS validation
				overlaps = cosmic_genome_tree.get_all_overlaps(curr_obj)

				for o_df in overlaps:
					cds = o_df['Mutation CDS']
					if nucSub in cds:
						AA_sub = o_df["Mutation AA"]
						AA_sub = AA_sub.replace("p.", "")
						newValues.append(AA_sub)

		newDict.update({cell : newValues})

	return newDict



def are_hits_in_cosmic(queryList, SNP_bool):
	""" given a set of genome position strings, searches for the ones
		that are in the COSMIC database. now supports SNPs and indels!! """
	ret = []

	if len(queryList) > 0:
		if SNP_bool: # exact string match
			curr_df = pd.DataFrame(queryList, columns=['posStr', 'ref', 'alt'])
			overlap = set(curr_df['posStr']).intersection(set(genomePos_laud_db))
			keep = curr_df.where(curr_df['posStr'].isin(overlap))
			keep = keep.dropna()

			ret = keep.values.tolist() # convert entire df to list

		else: # indel case
			ret = []
			for i in range(0,len(queryList)):
				pos_str = queryList[i]
				pos_str_raw = pos_str[0]
				curr_obj = utils.GenomePosition.from_str(pos_str_raw)
				try:
					b = cosmic_genome_tree.get_best_overlap(curr_obj)
				
					if b is not None:
						ret.append(pos_str)
				except AttributeError:
					continue

	return(ret)



def build_genome_positions_dict(fileName):
	""" creates dict with genome coords for cosmic filtered hits to specific GOI """

	cell = fileName.replace(cwd + "scVCF_filtered_all/", "")
	cell = cell.replace(".vcf", "")	

	df = vcf_to_dataframe(fileName)
	genomePos_query = df.apply(generate_genome_pos_str, axis=1)
	genomePos_query_breakdown = breakdown_by_mutation_type(genomePos_query)
	
	SNPs = genomePos_query_breakdown[0]
	indels = genomePos_query_breakdown[1]

	# get the entries shared between curr cells VCF and the cosmic filter set
	shared_SNPs = are_hits_in_cosmic(SNPs, 1)
	shared_indels = are_hits_in_cosmic(indels, 0)

	all_shared_regions = shared_SNPs + shared_indels

	ret = [cell, all_shared_regions]
	return(ret)



""" get cmdline input """
@click.command()
@click.option('--gene', default = 'EGFR', prompt='gene_id (all caps)', required=True, type=str)
@click.option('--nthread', default = 2, prompt='number of threads', required=True, type=int)
@click.option('--outprefix', default = 'sampleOut', prompt='prefix to use for outfile', required=True, type=str)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)
 


def get_aa_mutations(gene, nthread, outprefix, wrkdir):
	""" for a specific gene of interest, get the complete set of amino acid level mutations
		for each cell in dataset """
	global database
	global database_laud
	global genomePos_laud_db
	global cosmic_genome_tree
	global cwd

	cwd = wrkdir
	gene_name = gene

	print(' ')
	print('setting up COSMIC database...')
	database = pd.read_csv(cwd + "CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = get_laud_db(gene_name)
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	# init interval tree
	cosmic_genome_tree = utils.GenomeDataframeTree(lambda row: utils.GenomePosition.from_str(str(row["Mutation genome position"])), database_laud)

	fNames = get_filenames()

	print('searching for relevant vcf hits')
	p = mp.Pool(processes=nthread)
		
	try:  
		goiList = list(tqdm(p.imap(build_genome_positions_dict, fNames), total=len(fNames)))
	finally:
		p.close()
		p.join()

	cells_dict_GOI_coords = {} # convert to dict

	for item in goiList:
		cell = item[0]
		muts = item[1]
		
		toAdd = {cell:muts}
		cells_dict_GOI_coords.update(toAdd)

	goiDict_AA = get_corresponding_aa_subs(cells_dict_GOI_coords)

	write_csv(cells_dict_GOI_coords, cwd + outprefix + '.csv')
	write_csv(goiDict_AA, cwd + outprefix + '_AA.csv')
