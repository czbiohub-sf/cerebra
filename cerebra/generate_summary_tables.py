""" creates the final validaion table(s)
	by SAMPLE and by CELL """

import summarize_module
import pandas as pd
import numpy as np


""" get cmdline input """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def generate_summary_tables(test, wrkdir):
	""" generate by cell and by sample summary tables for your experiment """
	test_bool = test
	cwd = wrkdir

	# read in all of these by-gene amino-acid level mutation counts objs
	egfrPATH = cwd + 'egfr_germline_out_AA.csv'
	brafPATH = cwd + 'braf_germline_out_AA.csv'
	krasPATH = cwd + 'kras_germline_out_AA.csv'

	egfr_df = pd.read_csv(egfrPATH, header=None, names=['cell', 'mutations'])
	braf_df = pd.read_csv(brafPATH, header=None, names=['cell', 'mutations'])
	kras_df = pd.read_csv(krasPATH, header=None, names=['cell', 'mutations'])

	# first step is to generate the mutationsDF
	mutationsDF = pd.DataFrame(columns=['cell', 'brafMut', 'egfrMut', 'krasMut'])
	mutationsDF['cell'] = egfr_df['cell']
	mutationsDF['egfrMut'] = egfr_df['mutations'] # fill in EGFR first
	summarizeModule.mutations_df_fill_in('braf', braf_df, mutationsDF) 
	summarizeModule.mutations_df_fill_in('kras', kras_df, mutationsDF)

	# converting lists to strs. makes downstream analysis easier
	summarizeModule.remove_extra_characters_mutations_df('egfr', mutationsDF)
	summarizeModule.remove_extra_characters_mutations_df('braf', mutationsDF)
	summarizeModule.remove_extra_characters_mutations_df('kras', mutationsDF)

	# read in patientMetadata
	metaPATH = cwd + 'metadata_all_cells_4.10.19.csv'
	patientMetadata = pd.read_csv(metaPATH)
	patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd

	# init the summary table
	cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'mutations_found_EGFR', 'mutations_found_BRAF', 'mutations_found_KRAS', 'fusions_found', 'tumorCell_bool']
	summaryTable = pd.DataFrame(columns=cols)
	summaryTable['cell'] = mutationsDF['cell']

	# fill in various metadata cols
	summarizeModule.generic_summary_table_fill_in('patient_id', 'patient', summaryTable, patientMetadata)
	summarizeModule.generic_summary_table_fill_in('driver_gene', 'clinical_driver_gene', summaryTable, patientMetadata)
	summarizeModule.generic_summary_table_fill_in('driver_mutation', 'clinical_mutation', summaryTable, patientMetadata)

	# fill in mutations found col 
	summaryTable['mutations_found_EGFR'] = mutationsDF['egfrMut']
	summaryTable['mutations_found_KRAS'] = mutationsDF['krasMut']
	summaryTable['mutations_found_BRAF'] = mutationsDF['brafMut']

	# read in fusions dataframe, then fill in summaryTable
	fusionsDF = pd.read_csv(cwd + 'fusion_dataframe.csv')
	summarizeModule.fusions_fill_in(fusionsDF, summaryTable)

	# set up a col to translate 'raw' mutation calls to 'clinical'
	summaryTable['mutations_found_translated'] = ""
	summarizeModule.translated_muts_fill_in_egfr(summaryTable)
	summarizeModule.translated_muts_fill_in('KRAS', summaryTable)
	summarizeModule.translated_muts_fill_in('BRAF', summaryTable)
	summarizeModule.translated_muts_fill_in_fusions(summaryTable)

	# convert lists to strs, so i can get set -- probably not necessary
	summarizeModule.convert_to_string(summaryTable)

	# fill in clin_mut_found_bool 
	summarizeModule.clin_mut_found_fill_in(summaryTable)
	summarizeModule.clin_mut_found_fill_in_fus(summaryTable)

	# fill in tumorCellBool
	summarizeModule.tumor_cell_bool_fill_in(summaryTable, cwd)

	# get per-cell ROI coverage dfs
	braf_V600E_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('braf', 'V600E', cwd)
	egfr_L858R_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'L858R', cwd)
	egfr_exon19del_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'exon19del', cwd)
	egfr_exon20ins_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'exon20ins', cwd)
	egfr_G719X_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'G719X', cwd)
	egfr_L861Q_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'L861Q', cwd)
	egfr_S768I_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'S768I', cwd)
	egfr_T790M_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('egfr', 'T790M', cwd)
	kras_G12C_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('kras', 'G12C', cwd)
	kras_G13X_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('kras', 'G13X', cwd)
	kras_Q61X_cov_nonZero = summarizeModule.get_non_zero_cov_ROI('kras', 'Q61X', cwd)

	# fix up some of the wierd ones
	kras_G13X_cov_nonZero['depth_gvcf'][4202] = 34
	kras_Q61X_cov_nonZero['depth_gvcf'][6431] = 92
	egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip('[')
	egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip(']')
	egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip("'")

	# fill in ROI coverage
	summarizeModule.ROI_coverage_fill_in(braf_V600E_cov_nonZero, 'BRAF', 'V600E', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_G719X_cov_nonZero, 'EGFR', 'G719X', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_L858R_cov_nonZero, 'EGFR', 'L858R', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_L861Q_cov_nonZero, 'EGFR', 'L861Q', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_S768I_cov_nonZero, 'EGFR', 'S768I', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_T790M_cov_nonZero, 'EGFR', 'T790M', summaryTable)
	summarizeModule.ROI_coverage_fill_in(kras_G12C_cov_nonZero, 'KRAS', 'G12C', summaryTable)
	summarizeModule.ROI_coverage_fill_in(kras_G13X_cov_nonZero, 'KRAS', 'G13X', summaryTable)
	summarizeModule.ROI_coverage_fill_in(kras_Q61X_cov_nonZero, 'KRAS', 'Q61X', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_exon19del_cov_nonZero, 'EGFR', 'del19', summaryTable)
	summarizeModule.ROI_coverage_fill_in(egfr_exon20ins_cov_nonZero, 'EGFR', 'ins20', summaryTable)

	# trim it down
	summaryTable_trimmed = summaryTable[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'tumorCell_bool', 'mutations_found_translated']]
	summaryTable_trimmed.columns = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool', 'mutations_found']
	summaryTable_trimmed = summaryTable_trimmed[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'mutations_found', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool']]
	summaryTable_trimmed

	# add sample name col to summary table
	summaryTable_trimmed['sample_name'] = ''
	summarizeModule.generic_summary_table_fill_in('sample_name', 'sample_name', summaryTable_trimmed, patientMetadata)

	# write
	summaryTable_trimmed.to_csv(cwd + 'validationTable_cells.csv', index=False)

	# get min set of sample names
	relevantSamplesSet = set(summaryTable_trimmed['sample_name'])
	relevantSamplesList = list(relevantSamplesSet)
	relevantSamplesSeries = pd.Series(relevantSamplesList)

	# init validationTable_samples
	cols = ['sample', 'patient', 'driver_gene', 'driver_mutation', 'mutations_found', 'numCells', 'numTumorCells', 'numTumorCells_w_coverage_to_ROI', 'numTumorCells_clinMut_found']
	validationTable_samples = pd.DataFrame(columns=cols)
	validationTable_samples['sample'] = relevantSamplesSeries

	# fill in metadata fields
	summarizeModule.validation_table_metadata_fill_in('patient_id', 'patient', validationTable_samples, patientMetadata)
	summarizeModule.validation_table_metadata_fill_in('driver_gene', 'driver_gene', validationTable_samples, patientMetadata)
	summarizeModule.validation_table_metadata_fill_in('driver_mutation', 'driver_mutation', validationTable_samples, patientMetadata)

	# fill in mutations found
	muts_dict = summarizeModule.validation_table_dict_muts(validationTable_samples, summaryTable_trimmed)
	validationTable_samples['mutations_found'] = muts_dict.values()

	# fill in numTumorCells (various)
	tc_dict = summarizeModule.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'tumorCell_bool')
	tc_cov_dict = summarizeModule.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'coverage_to_ROI')
	clinMut_dict = summarizeModule.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'clinical_mutation_found_bool')

	validationTable_samples['numTumorCells'] = tc_dict.values()
	validationTable_samples['numTumorCells_w_coverage_to_ROI'] = tc_cov_dict.values()
	validationTable_samples['numTumorCells_clinMut_found'] = clinMut_dict.values()

	# clean up 
	validationTable_samples = validationTable_samples.drop([0]) # this can change
	cols = ['sample', 'patient', 'driver_gene', 'driver_mutation', 'mutations_found', 'numTumorCells', 'numTumorCells_w_coverage_to_ROI', 'numTumorCells_clinMut_found']
	validationTable_samples = validationTable_samples[cols]

	# write
	validationTable_samples.to_csv(cwd + 'validationTable_samples.csv', index=False)

