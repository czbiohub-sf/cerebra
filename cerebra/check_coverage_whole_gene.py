""" this module reproduces figure 2E """

import os
import gzip
import pandas as pd
import warnings
import click 
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



def GOI_df_subset(vcf_, chrom_, start_, end_):
	""" subset a single vcf based on genomic coords"""
	chrStr = 'chr' + str(chrom_)
    
	keep0 = vcf_['CHROM'] == chrStr
	vcf_sub0 = vcf_[keep0]

	keep1 = vcf_sub0['POS'] >= start_
	vcf_sub1 = vcf_sub0[keep1]

	keep2 = vcf_sub1['POS'] <= end_
	vcf_sub2 = vcf_sub1[keep2]
    
	return(vcf_sub2)



def coverage_search(df):
	""" given subset of vcf entries, search for AD (DepthPerAlleleBySample) col """ 
	counts = []
	for i in range(0, len(df.index)):
		row = df.iloc[i]
		extra_col = str(row['20'])

		try:
			AD = extra_col.split(':')[1]
			wt_count = int(AD.split(',')[0])
			variant_count = int(AD.split(',')[1])
			total_count = wt_count + variant_count

			ratio = str(variant_count) + ':' + str(total_count)
			counts.append(ratio)

		except: # picking up a wierd edge case -- '20' field is malformed
			print(extra_col)
			continue
		
	return(counts)



def driver(cellFile):
	""" for a given cell, search for variants to all genes! """
	currCell = cellFile.strip(wrkdir_ + 'scVCF_filtered_all/')
	currCell = currCell.strip('.vcf')
	#print(currCell)

	vcf = vcf_to_dataframe(cellFile)
	GOI_counts_by_cell = []

	for gene in gene_names:
		keep = hg38.gene_id == gene
		hg38_GOI_subset = hg38[keep]
		hg38_GOI_subset = hg38_GOI_subset.reset_index(drop=True)

		start = min(list(hg38_GOI_subset.start))
		end = max(list(hg38_GOI_subset.end))
		seqname = hg38_GOI_subset.seqname.iloc[0]
		chrom = seqname.strip('chr')

		vcf_sub = GOI_df_subset(vcf, chrom, start, end)
    
		if not vcf_sub.empty:
			counts = coverage_search(vcf_sub)
		else:
			counts = ['0:0']

		GOI_counts_by_cell.append([gene, counts])

	ret = [currCell, GOI_counts_by_cell]
	return(ret)



def average_by_gene(quad_list_, ratios_df_):
	""" averages expression across found loci for each gene, and converts to 
		variant / total read ratio """

	for elm in quad_list_:
		cell_name = elm[0]
		genes_breakdown = elm[1]

		for sub_elm in genes_breakdown:
			gene = sub_elm[0]
			vals = sub_elm[1]
			v_vals_total = 0
			t_vals_total = 0
			
			for val in vals:
				v_val = int(val.split(':')[0])
				t_val = int(val.split(':')[1])
				v_vals_total += v_val
				t_vals_total += t_val

			if t_vals_total != 0:
				ratio = v_vals_total / t_vals_total
				ratio = round(ratio, 2) # round to two decimal places
			else:
				ratio = 0
			
			ratios_df_[gene][cell_name] = ratio

	return(ratios_df_)



""" get cmdline input """
@click.command()
@click.option('--wrkdir', default = '/home/ubuntu/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)
@click.option('--genes_list', default = 'genesList.csv', prompt='name of csv file with genes of interest to evaluate coverage for. should be in wrkdir', required=True, type=str)
@click.option('--nthread', default = 16, prompt='number of threads', required=True, type=int)



def check_coverage_whole_gene(wrkdir, genes_list, nthread):
	""" evaluate the coverage across every loci with a reported variant
		for each gene in a user defined list. reports on a per-gene 
		basis. """

	global gene_names
	global hg38
	global wrkdir_
	num_proc = 16

	wrkdir_ = wrkdir # cmdline inputs cant be declared globals? 

	vcf_list = os.listdir(wrkdir_ + 'scVCF_filtered_all/')
	cell_files_list = [wrkdir_ + 'scVCF_filtered_all/' + s for s in vcf_list] # list comprehension
	cells_list = [elm.strip(wrkdir_ + 'scVCF_filtered_all' + '.vcf') for elm in cell_files_list]

	GOI_df = pd.read_csv(wrkdir_ + genes_list, header=None, names=['gene'])
	gene_names = list(GOI_df.gene)

	print(' ')
	hg38 = pd.read_csv(wrkdir_ + 'hg38-plus.gtf', sep='\t', header=None, names=['seqname', 
			'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

	hg38['gene_id'] = 'NA'
	gene_ids = hg38.attribute.apply(lambda elm: elm.split(';')[0].split('"')[1])
	hg38.gene_id = gene_ids


	print('creating pool')
	p = mp.Pool(processes=num_proc)
	print('running...')

	try:
		quad_list = p.map(driver, cell_files_list, chunksize=1) # returns quadruply nested list
	finally:
		p.close()
		p.join()

	ratios_df = pd.DataFrame(columns=gene_names, index=cells_list)
	ratios_df[:] = 0
	ratios_df = ratios_df.astype(float)

	ratios_df_filled = average_by_gene(quad_list, ratios_df)
	ratios_df_filled.to_csv(wrkdir_ + 'ratios_df_test.csv')
