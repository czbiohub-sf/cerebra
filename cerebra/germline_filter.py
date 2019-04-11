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

	cwd = os.getcwd()
	filterDir = cwd + '/wrkdir/filteredOut/'
	filterDir_list = os.listdir(filterDir)

	filteredCells = []
	for f in filterDir_list:
		cell = f.strip('_unique.vcf')
		filteredCells.append(cell)

	epiDir = cwd + '/wrkdir/scVCF/'
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

	os.system('sudo mkdir -p wrkdir/scVCF_filtered_all/')
	os.system('sudo chmod -R 777 wrkdir/scVCF_filtered_all/')

	# copy over the non-filtered cells
	outPATH = cwd + '/wrkdir/scVCF_filtered_all/'
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
	with open('vcfheader.txt', 'r') as f:
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
	basePATH = os.getcwd()
	patientPATH = basePATH + '/wrkdir/bulkVCF/' + patient
	cellPATH = basePATH + '/wrkdir/scVCF/' + cell + '.vcf'
	
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



""" get cmdline input """
@click.command()
@click.option('--metadata', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/metadata_all_cells_4.10.19.csv', prompt='s3 path to by-cell metadata', required=True, type=str)
@click.option('--sc_vcf', default = 's3://lincoln.harris-work/scVCF/', prompt='s3 path to single-cell VCFs', required=True, type=str)
@click.option('--bulk_vcf', default = 's3://lincoln.harris-work/bulkVCF/', prompt='s3 path to bulk VCFs', required=True, type=str)
@click.option('--outpath', default = 's3://lincoln.harris-work/filteredOut/', prompt='s3 path to where filtered output files should be pushed', required=True, type=str)



def germline_filter(metadata, sc_vcf, bulk_vcf, outpath):
	""" driver function. """
	global patientMetadata

	# read in patient metadata
	#os.system('sudo mkdir -p wrkdir')
	#os.system('sudo chmod -R 777 wrkdir')
	#cmd = 'aws s3 cp ' + metadata + ' wrkdir/ --quiet'
	#os.system(cmd)
	patientMetadata = pd.read_csv('wrkdir/metadata_all_cells_4.10.19.csv')

	# get a list of all the single-cell VCF files
	#os.system('sudo mkdir -p wrkdir/scVCF')
	#os.system('sudo chmod -R 777 wrkdir/scVCF')
	#cmd = 'aws s3 cp ' + sc_vcf + ' wrkdir/scVCF/ --recursive --quiet'
	#os.system(cmd)
	cwd = os.getcwd()
	vcfDir = cwd + '/wrkdir/scVCF/'
	scVCF_list = os.listdir(vcfDir)

	# get list of bulk VCF files
	#os.system('sudo mkdir -p wrkdir/bulkVCF')
	#os.system('sudo chmod -R 777 wrkdir/bulkVCF/')
	#cmd = 'aws s3 cp ' + bulk_vcf + ' wrkdir/bulkVCF/ --recursive --quiet'
	#os.system(cmd)
	bulkVCF_dir = cwd + '/wrkdir/bulkVCF/'
	bulkVCF_list = os.listdir(bulkVCF_dir)

	patientsRun = []

	os.system('sudo mkdir -p wrkdir/filteredOut')
	os.system('sudo chmod -R 777 wrkdir/filteredOut/')

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
				outStr = 'wrkdir/filteredOut/' + currCell + '_unique.vcf'
				write_vcf(currCell_unique, outStr)
			
			patientsRun.append(currPatient)

	create_final_outdir()

	# push outfiles, then clean up
	#cmd = 'aws s3 cp wrkdir/filteredOut/ ' + outpath + ' --recursive --quiet'
	#os.system(cmd)
	#os.system('rm -rf wrkdir')
