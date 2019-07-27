""" evaluate coverage on a loci-specific basis. most of these functions are
	copied from get_aa_mutations.py """

import numpy as np
import os
from . import utils
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
	for file in os.listdir(cwd + "vcf_test_set/"):
		if file.endswith(".vcf"):
			fullPath = cwd + 'vcf_test_set/' + file 
			files.append(fullPath)
    
	return files



def get_laud_db(gene_, db):
    """ returns the COSMIC database after lung and fathmm filter  """
    pSiteList = db.index[db['Primary site'] == 'lung'].tolist()
    db_filter = db.iloc[pSiteList]
    #keepRows = db_filter['FATHMM score'] >= 0.7
    #db_fathmm_filter = dbfilter[keepRows]
    #db_fathmm_filter = db_fathmm_filter.reset_index(drop=True)
    keep = db_filter['Gene name'] == gene_
    db_gene = db_filter[keep]
    db_gene = db_gene.reset_index(drop=True)

    if len(db_gene.index) == 0:
    	print('   this gene is not in the cosmic database')
    	print('   	are you using the official HGNC gene name? ')
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



def get_rev_comp(b):
	""" return the reverse complement of a given base """
	if b == 'A':
		return('T')
	if b == 'T':
		return('A')
	if b == 'C':
		return('G')
	if b == 'G':
		return('C')



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
	cell = fileName.replace(cwd + "vcf_test_set/", "")
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



def GOI_df_subset(vcf_, chrom_, start_, end_):
	""" subset a single vcf based on genomic coords"""
	chrStr = 'chr' + str(chrom_)
	start_ = int(start_)
	end_ = int(end_)

	keep0 = vcf_['CHROM'] == chrStr
	vcf_sub0 = vcf_[keep0]

	keep1 = vcf_sub0['POS'] >= start_
	vcf_sub1 = vcf_sub0[keep1]

	keep2 = vcf_sub1['POS'] <= end_
	vcf_sub2 = vcf_sub1[keep2]
    
	return(vcf_sub2)



def coverage_search_on_vcf(df):
	""" given subset of vcf entries, search for AD (DepthPerAlleleBySample) col """ 
	counts = ''
	for i in range(0, len(df.index)):
		row = df.iloc[i]
		extra_col = str(row['20'])

		try:
			AD = extra_col.split(':')[1]
			wt_count = int(AD.split(',')[0])
			variant_count = int(AD.split(',')[1])
			total_count = wt_count + variant_count

			ratio = str(variant_count) + ':' + str(total_count)
			counts = ratio

		except: # picking up a wierd edge case -- '20' field is malformed
			print(extra_col)
			continue
		
	return(counts)



def get_corresponding_aa_sub(position_sub_str):
	""" for a given chromosomal position, searches the cosmic interval tree 
		for that exact position, and returns the corresponding AA level substitution """

	AA_sub = '?' # base case -- hits NOT in cosmic
	posStr = position_sub_str[0]
	ref = position_sub_str[1]
	alt = position_sub_str[2]

	curr_obj = utils.GenomePosition.from_str(posStr)

	if len(ref) > 1 or len(alt) > 1: # indel case -- just return best overlap
		b = cosmic_genome_tree.get_best_overlap(curr_obj)
		if b is not None:
			AA_sub = b["Mutation AA"]
			AA_sub = AA_sub.replace("p.", "")

	else: # SNP case -- specific CDS validation
		overlaps = cosmic_genome_tree.get_all_overlaps(curr_obj)
		nucSub = ref + '>' + alt

		for o_df in overlaps:
			cds = o_df['Mutation CDS']
			strand = o_df['Mutation strand']
			# taking into account cosmic's wierd strandedness field
			if strand == '-':
				nucSub = get_rev_comp(ref) + '>' + get_rev_comp(alt)

			if nucSub in cds:
				AA_sub = o_df["Mutation AA"]
				AA_sub = AA_sub.replace("p.", "")

	return(AA_sub)



def evaluate_coverage_driver(ROI_hits_dict, gene_, cd):
	""" takes a dict of ROI strings and converts to AA level muts, then
		calls coverage_search_on_vcf() for each of those AA level hitsß """

	for cell in ROI_hits_dict.keys():
		vcf_path = cwd + 'vcf_test_set/' + cell + '.vcf'
		vcf = vcf_to_dataframe(vcf_path)

		ROIs = ROI_hits_dict.get(cell)
		if len(ROIs) > 0:

		 	for hit in ROIs:
		 		aa_sub = get_corresponding_aa_sub(hit)
		 		posStr = hit[0]
		 		chrom = posStr.split(':')[0]
		 		start = posStr.split('-')[0].split(':')[1]
		 		end = posStr.split('-')[1]

		 		vcf_sub = GOI_df_subset(vcf, chrom, start, end)
		 		counts = coverage_search_on_vcf(vcf_sub)

		 		if cell in cd:
		 			curr_val = cd.get(cell)
		 			curr_val.append([gene_ + '_' + aa_sub, counts])
		 			cd.update({cell:curr_val})
		 		else:
		 			to_add = {cell:[[gene_ + '_' + aa_sub, counts]]}
		 			cd.update(to_add)

	return(cd)



def convert_to_df(cd):
	""" takes in a dictionary obj where keys are cells and values are loci and their
		associated coverage ratios, and converts to a cell x mutation dataframe """

	cells = list(cd.keys())
	l = list(cd.values())
	muts_list = []

	for item in l:
		for sub_item in item:
			curr_mut = sub_item[0]
			if '?' not in curr_mut:
				muts_list.append(curr_mut)

	muts_list = list(set(muts_list))
	df = pd.DataFrame(columns=muts_list, index=cells)
	df[:] = '0:0'

	for cell in cd.keys(): # write to coverage_dataframe
		values = cd.get(cell)
		for val in values:
			mut = val[0]
			cov = val[1]
			if '?' not in mut:
				df.loc[cell,mut] = cov

	# alphabetize rows and cols
	df = df.reindex(sorted(df.columns), axis=1)
	df = df.reindex(sorted(df.index), axis=0)
	df.index.name = 'sample'

	return(df)
 


""" get cmdline input """
@click.command()
@click.option('--genes_list', default = 'genesList.csv', prompt='name of csv file with genes of interest to evaluate coverage for. should be in wrkdir', required=True, type=str)
@click.option('--nthread', default = 2, prompt='number of threads', required=True, type=int)
@click.option('--outprefix', default = 'sampleOut.csv', prompt='prefix to use for outfile', required=True, type=str)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)
 


def check_coverage_loci(genes_list, nthread, outprefix, wrkdir):
	""" evaluate coverage for each loci for which we find a variant, for 
		a given set of genes  """
	global genomePos_laud_db
	global cosmic_genome_tree
	global cwd

	cwd = wrkdir
	fNames = get_filenames()
	GOI_df = pd.read_csv(cwd + genes_list, header=None, names=['gene'])
	gene_names = list(GOI_df.gene)
	database = pd.read_csv(cwd + "CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	coverage_dict = {}

	# driver loop 
	for gene in gene_names:
		print(gene)
		database_laud = get_laud_db(gene, database)
		genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

		# init interval tree
		cosmic_genome_tree = utils.GenomeIntervalTree(
            lambda row: utils.GenomePosition.from_str(str(row["Mutation genome position"])),
            (record for idx, record in database_laud.iterrows()))

		p = mp.Pool(processes=nthread)
			
		try:  
			goiList = list(p.imap(build_genome_positions_dict, fNames))
		finally:
			p.close()
			p.join()

		cells_dict_GOI_coords = {} # convert to dict

		for item in goiList:
			cell = item[0]
			muts = item[1]
			
			toAdd = {cell:muts}
			cells_dict_GOI_coords.update(toAdd)

		coverage_dict = evaluate_coverage_driver(cells_dict_GOI_coords, gene, coverage_dict)
	
	coverage_df = convert_to_df(coverage_dict)
	coverage_df.to_csv(cwd + outprefix)
