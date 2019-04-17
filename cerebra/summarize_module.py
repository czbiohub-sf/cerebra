""" module of functions that allow you to create per-cell / per-sample summary tables """

import numpy as np
import math
import pandas as pd


def mutations_df_fill_in(GOI, GOI_df, mutationsDF_):
	""" creates a cell-wise dataframe with mutations to each GOI """
	mutName = GOI + 'Mut'
	for i in range(0,len(mutationsDF_.index)):
		currCell = mutationsDF_['cell'][i]

		rightIndex = GOI_df['cell'] == currCell
		rightRow = GOI_df[rightIndex]
    
		rightCell = rightRow['cell']
		rightCell = str(rightCell).split()[1]
    
		rightMut = rightRow['mutations']
		rightMut = str(rightMut).split()[1]
    
		mutationsDF_[mutName][i] = rightMut



def remove_extra_characters_mutations_df(GOI, mutationsDF_):
	""" converting df cols from lists to strings """
	mutName = GOI + 'Mut'

	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("'", "") # remove quotes
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("[", "") # remove brackets
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("]", "") # remove brackets
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace(" ", "") # remove whitespace?



def generic_summary_table_fill_in(metaField, summaryField, summaryTable_, patientMetadata_):
	""" fills in a given metadata field in summaryTable_ """
	for i in range(0,len(summaryTable_.index)):
		currCell = summaryTable_['cell'].iloc[i]
		currPlate = currCell.split('_')[1]
    
		index_to_keep = patientMetadata_['plate'] == currPlate
		keepRow = patientMetadata_[index_to_keep]
		try:
			currField = list(keepRow[metaField])[0]
			summaryTable_[summaryField][i] = currField
		except IndexError:
			continue
			#print('ERROR: plate not found') # these are just the plates were NOT 
        	                                 # including in the analysis



def fusions_fill_in(fusionsDF_, summaryTable_):
	""" takes the existing fusionsDF and populates summaryTable_ with this shit """
	for i in range(0, len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		fusionsListCurr = []
    
		colList0 = list(fusionsDF_['ALK--EML4'])
		colList1 = list(fusionsDF_['ALK_any'])
		colList2 = list(fusionsDF_['EML4_any'])
		colList3 = list(fusionsDF_['NTRK_any'])
		colList4 = list(fusionsDF_['RET_any'])
		colList5 = list(fusionsDF_['ROS1_any'])

		if currCell in colList0:
			fusionsListCurr.append('ALK-EML4')
		elif currCell in colList1:
			fusionsListCurr.append('ALK_any')
		elif currCell in colList2:
			fusionsListCurr.append('EML4_any')
		elif currCell in colList3:
			fusionsListCurr.append('NTRK_any')
		elif currCell in colList4:
			fusionsListCurr.append('RET_any')
		elif currCell in colList5:
			fusionsListCurr.append('ROS1_any')
		else:
			fusionsListCurr = ""
        
		fusionsListCurr = str(fusionsListCurr)
		fusionsListCurr = fusionsListCurr.strip(']')
		fusionsListCurr = fusionsListCurr.strip('[')
		fusionsListCurr = fusionsListCurr.strip("'")
		fusionsListCurr = fusionsListCurr.strip(" ")
 
		summaryTable_['fusions_found'][i] = fusionsListCurr


def translated_muts_fill_in(GOI, summaryTable_):
	""" converts 'raw' mutation calls to something that more resembles
		those reported in our clinical cols. general """
	colName = 'mutations_found_' + GOI
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currMuts = summaryTable_[colName].iloc[i]
		currMuts_split = currMuts.split(',')
		for item in currMuts_split:
			if item != '' and '?' not in item:
				translatedList.append(GOI + ' ' + item)

		summaryTable_['mutations_found_translated'][i] = summaryTable_['mutations_found_translated'][i] + translatedList



def translated_muts_fill_in_egfr(summaryTable_):
	""" converts 'raw' mutation calls to something that more resembles
		those reported in our clinical cols. egfr, specificially """
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currMuts_egfr = summaryTable_['mutations_found_EGFR'].iloc[i]
		currMuts_egfr_split = currMuts_egfr.split(',')
		for item in currMuts_egfr_split:
			if 'delELR' in item:
				translatedList.append('EGFR del19')
			elif '745_' in item:
				translatedList.append('EGFR del19')
			elif '746_' in item:
				translatedList.append('EGFR del19')
			elif 'ins' in item:
				translatedList.append('EGFR ins20')
			elif item != '':
				translatedList.append('EGFR ' + item)
        
		summaryTable_['mutations_found_translated'][i] = translatedList



def translated_muts_fill_in_fusions(summaryTable_):
	""" converts 'raw' mutation calls to something that more resembles
		those reported in our clinical cols. for fusions """
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currFus = summaryTable_['fusions_found'].iloc[i]
		currFus_split = currFus.split(',')
		for item in currFus_split:
			if item == 'ALK-EML4':
				translatedList.append('ALK fusion')
				translatedList.append('EML4 fusion')
				translatedList.append('ALK-EML4 fusion')
			elif item != '' and '?' not in item:
				item = item.split('_')[0]
				translatedList.append(item + ' fusion')

		summaryTable_['mutations_found_translated'][i] = summaryTable_['mutations_found_translated'][i] + translatedList



def convert_to_string(summaryTable_):
	""" converting mutations_found_translated col from list to str. """
	for i in range(0,len(summaryTable_.index)):
		currStr = str(summaryTable_['mutations_found_translated'][i])
		currStr = currStr.replace("'", "")
		currStr = currStr.replace("]", "")
		currStr = currStr.replace("[", "")
		summaryTable_['mutations_found_translated'][i] = currStr



def clin_mut_found_fill_in(summaryTable_):
	""" fills in clin_mut_found_bool col: 1 if clin mut found, 0 if else """
	for i in range(0,len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		currMuts = summaryTable_['mutations_found_translated'][i]
		currClinGene = summaryTable_['clinical_driver_gene'][i]
		currClinMut = summaryTable_['clinical_mutation'][i]
		currClinMut_str = str(currClinGene) + ' ' + str(currClinMut)
    
		if currClinMut_str in currMuts:
			summaryTable_['clin_mut_found_bool'][i] = 1
		else:
			summaryTable_['clin_mut_found_bool'][i] = 0



def clin_mut_found_fill_in_fus(summaryTable_):
	""" fills in clin_mut_found_bool col: 1 if clin mut found, 0 if else
		but for fusions """
	for i in range(0,len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		currFus = summaryTable_['fusions_found'][i]
		currFus = currFus.strip('_any')
		currFus = currFus.split('-')[0]

		summaryTable_['clin_mut_found_bool'][i] = 0
		currClinGene = summaryTable_['clinical_driver_gene'][i]

		if currClinGene == currFus:
			summaryTable_['clin_mut_found_bool'][i] = 1



def tumor_cell_bool_fill_in(summaryTable_, cwd_):
	""" 1 if were calling the cell TUMOR in our seurat obj, 
		0 if else """
	# read in Seurat metadata
	metaPATH = cwd_ + 'metadataSeurat.csv'
	metadataSeurat = pd.read_csv(metaPATH)

	myCols = list(metadataSeurat.columns)
	myCols[0] = 'cell'
	metadataSeurat.columns = myCols
	
	indicies = metadataSeurat['inferCNV_annotation'] == 'perturbed'
	metadataSeurat_pert = metadataSeurat[indicies]
	
	tumorCellsList = list(metadataSeurat_pert['cell'])

	# now fill in 'tumorCell_bool' for summaryTable_
	for i in range(0, len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		if currCell in tumorCellsList:
			summaryTable_['tumorCell_bool'][i] = 1
		else:
			summaryTable_['tumorCell_bool'][i] = 0



def get_non_zero_cov_ROI(gene, mut, cwd_): 
	""" removes non-zero vals from given coverageByCell dataframe """
	fPATH = cwd_ + 'coverage/' + gene + '_' + mut + '_' + 'coverageByCell.csv'
	cov = pd.read_csv(fPATH)
	indices = cov['depth_gvcf'] != 0
	cov_nonZero = cov[indices]

	return(cov_nonZero)



def ROI_coverage_fill_in(coverage_df, queryGene, queryMutation, summaryTable_):
	""" fills in coverage for a given ROI, for summaryTable_ """
	for i in range(0, len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		currDriver = summaryTable_['clinical_driver_gene'][i]
		currMut = summaryTable_['clinical_mutation'][i]
    
		if currDriver == queryGene and currMut == queryMutation:
			if currCell in list(coverage_df['cellName']):
				index_cov_nonZero = coverage_df['cellName'] == currCell
				currRow_cov_nonZero = coverage_df[index_cov_nonZero]
				currDepth_gvcf = int(currRow_cov_nonZero['depth_gvcf'])
        
				summaryTable_['coverage_to_ROI'][i] = currDepth_gvcf
			else:
				summaryTable_['coverage_to_ROI'][i] = 0


   
def validation_table_metadata_fill_in(metaField, validationField, validationTable_, patientMetadata_):
	""" fills in metadata field for validationTable_ """
	for i in range(0, len(validationTable_.index)):
		currSample = validationTable_['sample'][i]
		try:
			rowToKeep = patientMetadata_['sample_name'] == currSample
			patientRows = patientMetadata_[rowToKeep] # will return MULTIPLE rows
			patientRows = patientRows.reset_index(drop=True)

			fillField = patientRows[metaField][0]
       
			validationTable_[validationField][i] = fillField
		except:
			continue
			#print('ERROR')



def validation_table_dict_muts(validationTable_, summaryTable_):
	""" returns dict that holds vals for all the muts to a 
		given cell """
	d = {}
	samplesList = validationTable_['sample']

	for item in samplesList:
		d.update({item:''})

	for i in range(0, len(summaryTable_.index)):
		currSample = summaryTable_['sample_name'][i]
		currMuts = summaryTable_['mutations_found'][i]
		currMuts = str(currMuts)
		currMutsSplit = currMuts.split(',')

		currDictVal = d[currSample]
    
		for item in currMutsSplit:
			if item not in currDictVal and item != 'nan':
				updateVal = currDictVal + item + ', '
				d.update({currSample:updateVal})

	return(d)


  
def validation_table_dict_generic(validationTable_, summaryTable_, field):
	""" returns dict that holds values for num cells that are tumor OR
		have coverage to a given ROI """
	d = {}
	samplesList = validationTable_['sample']
	for item in samplesList:
		d.update({item:0})

	for i in range(0, len(summaryTable_.index)):
		currSample = summaryTable_['sample_name'][i]
		currBool = summaryTable_[field][i]

		currDictVal = d[currSample]  

		if not math.isnan(currBool) and currBool != 0:
			updateVal = currDictVal + 1
			d.update({currSample:updateVal})

	return(d)
    