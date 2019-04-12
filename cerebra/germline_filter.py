"""
This program takes in two sets of vcfs, single-cell and bulk (peripheral
blood, ie. germline) and filters out the common variants found in both sc 
and bulk. Creates a new directory, filteredOut, that contans the filtered vcfs. 
"""

import os
import sys
import shutil
import warnings
import pandas as pd
import click
import VCF
warnings.simplefilter(action='ignore', category=FutureWarning)


def create_final_outdir():
	""" creates a finalized out dir that has all of the filtered vcfs as
		well as the ones we didnt have germline for """

	filterDir = cwd + 'filteredOut/'
	filterDir_list = os.listdir(filterDir)

	filteredCells = []
	for f in filterDir_list:
		cell = f.strip('_unique.vcf')
		filteredCells.append(cell)

	epiDir = cwd + 'scVCF/'
	epiDir_list = os.listdir(epiDir)

	epiCells = []
	for f in epiDir_list:
		cell = f.strip('.vcf')
		epiCells.append(cell)

	# get cells in epiCells but NOT filteredCells
	nonFilteredCells = set(epiCells) - set(filteredCells)

	nonFilteredCells_list = []
	for cell in nonFilteredCells:
		f = cell + '.vcf' 
		nonFilteredCells_list.append(f)

	cmd1 = 'sudo mkdir -p ' + cwd + 'scVCF_filtered_all/'
	cmd2 = 'sudo chmod -R 777 ' + cwd + 'scVCF_filtered_all/'
	os.system(cmd1)
	os.system(cmd2)

	# copy over the non-filtered cells
	outPATH = cwd + 'scVCF_filtered_all/'
	for file in nonFilteredCells_list:
		src = epiDir + file
		dst = outPATH + file
		shutil.copyfile(src, dst)

	# copy over all the filtered cells
	for file in filterDir_list:
		f = file.strip('_unique.vcf')
		f = f + '.vcf'
		src = filterDir + file
		dst = outPATH + f
		shutil.copyfile(src, dst)



def write_vcf(df, outStr_):
	""" routine for writing VCF files, from an existing dataframe. 
		essentially just adding in this horrible vcf header. """
	with open(cwd + 'vcfheader.txt', 'r') as f:
		header = f.read()
	
		df['QUAL'] = 60.28
		df['FILTER'] = '.'
		df['INFO'] = 'AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=3.00;QD=30.14;SOR=2.303'

		output_VCF = outStr_
		with open(output_VCF, 'w') as vcf:
			vcf.write(header)

		df.to_csv(output_VCF, sep="\t", mode='a', index=False)



def get_patient_cells_list(scVCF_list_, patientID):
	""" get the list of cell names from a given patient """
	currPatient_cells_ = []

	for item in scVCF_list_:
		currCell = item.strip('.vcf')
		currPlate = currCell.split('_')[1]
		rowToKeep = patientMetadata['plate'] == currPlate
    
		try:
			currPatient = patientMetadata['patient_id'][rowToKeep]
			index = currPatient.index[0]
			currPatientVal = currPatient[index]
			if currPatientVal == patientID:
				currPatient_cells_.append(currCell)
		except:
			continue
	print('numCells: %d' % len(currPatient_cells_))
	return currPatient_cells_



def get_unique_vcf_entries(patient, cell):
	""" do the germline filter, and return a dataframe with only the 
		UNIQUE entries for a given cell """
	patientPATH = cwd + 'bulkVCF/' + patient
	cellPATH = cwd + 'scVCF/' + cell + '.vcf'
	
	try:
		patient_df = VCF.dataframe(patientPATH)
		cell_df = VCF.dataframe(cellPATH)
	except FileNotFoundError:
		print('FILE NOT FOUND: %s' % cellPATH)
		return
    
	patient_df_trimmed = patient_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
	cell_df_trimmed = cell_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
    
	# get whats SHARED between patient and cell 
	patient_cell_concat = pd.concat([patient_df_trimmed, cell_df_trimmed])
	rowsToKeep = patient_cell_concat.duplicated()
	patient_cell_shared = patient_cell_concat[rowsToKeep]
	patient_cell_shared = patient_cell_shared.reset_index(drop=True)

	# now go back to the original cell df, pull out whats UNIQUE 
	cell_cell_concat = pd.concat([cell_df_trimmed, patient_cell_shared])
	cell_cell_concat_noDups = cell_cell_concat.drop_duplicates(keep=False)
	cell_cell_concat_noDups = cell_cell_concat_noDups.reset_index(drop=True)
    
	return(cell_cell_concat_noDups)



""" launch """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def germline_filter(test, wrkdir):
	""" given a set of single-cell vcfs and bulk-seq vcfs (peripheral blood), this
		program subtracts out the mutations common to sc- and bulkVCF. """
	global patientMetadata
	global cwd

	cwd = wrkdir

	# read in patient metadata
	patientMetadata = pd.read_csv(cwd + 'metadata_all_cells_4.10.19.csv')

	# get a list of all the single-cell VCF files
	vcfDir = cwd + 'scVCF/'
	scVCF_list = os.listdir(vcfDir)

	# get list of bulk VCF files
	bulkVCF_dir = cwd + 'bulkVCF/'
	bulkVCF_list = os.listdir(bulkVCF_dir)

	patientsRun = []

	cmd1 = 'sudo mkdir -p ' + cwd + 'filteredOut'
	cmd2 = 'sudo chmod -R 777 ' + cwd + 'filteredOut/' 
	os.system(cmd1)
	os.system(cmd2)

	# outer loop -- by PATIENT
	for item in bulkVCF_list:
		currSample = item.strip('.vcf')
		currPatient = currSample.split('_')[0]
		suffix1 = currSample.split('_')[1]
		try:
			suffix2 = currSample.split('_')[2]
		except IndexError:
			suffix2 = ''
	
		if suffix2 != '' and currPatient not in patientsRun:
			print('WHOLE BLOOD FOUND, for %s' % currPatient)
			currPatient_cells = get_patient_cells_list(scVCF_list, currPatient)

			# inner loop -- by CELL 
			for currCell in currPatient_cells:
				currCell_unique = get_unique_vcf_entries(item, currCell)
				outStr = cwd + 'filteredOut/' + currCell + '_unique.vcf'
				write_vcf(currCell_unique, outStr)
			
			patientsRun.append(currPatient)

	if test:
		cmd1 = 'sudo mkdir ' + cwd + 'test_germline_filter'
		cmd2 = 'sudo chmod -R 777 ' + cwd + 'test_germline_filter'
		cmd3 = 'cp -r ' + cwd + 'filteredOut/ ' + cwd + 'test_germline_filter/' 
		os.system(cmd1)
		os.system(cmd2)
		os.system(cmd3)
	else:
		create_final_outdir()
