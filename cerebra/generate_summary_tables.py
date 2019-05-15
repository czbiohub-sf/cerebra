""" creates the final validaion table(s)
	by SAMPLE and by CELL """

import summarize_module
import pandas as pd
import numpy as np
import click


""" get cmdline input """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/home/ubuntu/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def generate_summary_tables(test, wrkdir):
	""" generate by cell and by sample summary tables for your experiment """
	test_bool = test
	cwd = wrkdir
	muts_path = cwd + 'geneSearch_out/'

	mutationsDF = pd.DataFrame(columns=['cell', 'AKT1_mut', 'ALK_mut', 'BAP1_mut', 'BRAF_mut', 'DDR2_mut', 
		'DROSHA_mut', 'EGFR_mut', 'ERBB2_mut', 'ERBB4_mut', 'FGFR2_mut', 'GRIN2A_mut', 'HIF1a_mut', 'KDR_mut', 
        'KEAP1_mut', 'KRAS_mut', 'MAP2K1_mut', 'MAP2K2_mut', 'MYCL_mut', 'NFE2L2_mut', 'NKX21_mut', 'NOTCH1_mut',
        'PIK3CB_mut', 'PTPN13_mut', 'PTPRT_mut', 'RAD21_mut', 'RB1_mut', 'RBM10_mut', 'SMARCA4_mut', 
        'SOX2_mut', 'STK11_mut', 'TP53_mut', 'TP63_mut', 'BRIP_mut', 'NF1_mut', 'BRCA2_mut', 'VHL_mut', 'BTK_mut', 
        'RARA_mut', 'SPTA1_mut', 'FGFR1_mut', 'SETD2_mut', 'JAK1_mut', 'SMAD4_mut', 'CTNNB1_mut', 'TET2_mut', 
        'PIK3CA_mut'])

	genesList = ['BRIP', 'AKT1', 'ALK', 'BAP1', 'BRAF', 'DDR2', 'DROSHA', 'EGFR', 'ERBB2', 'ERBB4', 'FGFR2', 'GRIN2A',
				'HIF1a', 'KDR', 'KEAP1', 'KRAS', 'MAP2K1', 'MAP2K2', 'MYCL', 'NFE2L2', 'NKX21', 'NOTCH1', 'PIK3CB',
				'PTPN13', 'PTPRT', 'RAD21', 'RB1', 'RBM10', 'SMARCA4', 'SOX2', 'STK11', 'TP53', 'TP63','NF1', 
				'BRCA2', 'VHL', 'BTK', 'RARA', 'SPTA1', 'FGFR1', 'SETD2', 'JAK1', 'SMAD4', 'CTNNB1', 'TET2', 'PIK3CA']

	# fill in EGFR first
	EGFR_path = muts_path + 'EGFR_AA.csv'
	EGFR_df = pd.read_csv(EGFR_path, header=None, names=['cell', 'mutations'])
	mutationsDF['cell'] = EGFR_df['cell']
	mutationsDF['EGFR_mut'] = EGFR_df['mutations'] 

	# now fill in everything else
	for gene in genesList:
		gene_path = muts_path + gene + '_AA.csv'
		gene_df = pd.read_csv(gene_path, header=None, names=['cell', 'mutations'])
		summarize_module.mutations_df_fill_in(gene, gene_df, mutationsDF)
		summarize_module.remove_extra_characters_mutations_df(gene, mutationsDF)

	# read in patientMetadata
	metaPATH = cwd + 'metadata_all_cells_4.10.19.csv'
	patientMetadata = pd.read_csv(metaPATH)
	patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd

	# init the summary table
	cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 
        'mutations_found_AKT1', 'mutations_found_ALK', 'mutations_found_BAP1', 'mutations_found_BRAF', 
        'mutations_found_DDR2', 'mutations_found_DROSHA', 'mutations_found_EGFR', 'mutations_found_ERBB2',
        'mutations_found_ERBB4', 'mutations_found_FGFR2', 'mutations_found_GRIN2A','mutations_found_HIF1a', 
        'mutations_found_KDR', 'mutations_found_KEAP1', 'mutations_found_KRAS', 'mutations_found_MAP2K1', 
        'mutations_found_MAP2K2', 'mutations_found_MYCL', 'mutations_found_NFE2L2', 'mutations_found_NKX21',
        'mutations_found_NOTCH1','mutations_found_PIK3CB', 'mutations_found_PTPN13', 'mutations_found_PTPRT',
        'mutations_found_RAD21', 'mutations_found_RB1', 'mutations_found_RBM10',
        'mutations_found_SMARCA4', 'mutations_found_SOX2', 'mutations_found_STK11', 'mutations_found_TP53',
        'mutations_found_TP63', 'mutations_found_BRIP', 'mutations_found_NF1', 'mutations_found_BRCA2', 
        'mutations_found_VHL', 'mutations_found_BTK', 'mutations_found_RARA', 'mutations_found_SPTA1', 
        'mutations_found_FGFR1', 'mutations_found_SETD2', 'mutations_found_JAK1', 'mutations_found_SMAD4',
        'mutations_found_CTNNB1', 'mutations_found_TET2', 'mutations_found_PIK3CA', 'fusions_found', 'tumorCell_bool']

	summaryTable = pd.DataFrame(columns=cols)
	summaryTable['cell'] = mutationsDF['cell']

	# fill in various metadata cols
	summarize_module.generic_summary_table_fill_in('patient_id', 'patient', summaryTable, patientMetadata)
	summarize_module.generic_summary_table_fill_in('driver_gene', 'clinical_driver_gene', summaryTable, patientMetadata)
	summarize_module.generic_summary_table_fill_in('driver_mutation', 'clinical_mutation', summaryTable, patientMetadata)

	# fill in mutations found col 
	for gene in genesList:
		summary_col = 'mutations_found_' + gene
		mutsdf_col = gene + '_mut'
		summaryTable[summary_col] = mutationsDF[mutsdf_col]

	# read in fusions dataframe, then fill in summaryTable
	fusionsDF = pd.read_csv(cwd + 'fusion_dataframe_4.19.csv')
	summarize_module.fusions_fill_in(fusionsDF, summaryTable)

	# set up a col to translate 'raw' mutation calls to 'clinical'
	summaryTable['mutations_found_translated'] = ""
	summarize_module.translated_muts_fill_in_egfr(summaryTable)

	# translate the rest of 'em
	for gene in genesList:
		summarize_module.translated_muts_fill_in(gene, summaryTable)

    # and the fusions
	summarize_module.translated_muts_fill_in_fusions(summaryTable)

	# convert lists to strs, so i can get set -- probably not necessary
	summarize_module.convert_to_string(summaryTable)

	# fill in clin_mut_found_bool 
	summarize_module.clin_mut_found_fill_in(summaryTable)
	summarize_module.clin_mut_found_fill_in_fus(summaryTable)

	# fill in tumorCellBool
	summarize_module.tumor_cell_bool_fill_in(summaryTable, cwd)

	# get per-cell ROI coverage dfs
	#braf_V600E_cov_nonZero = summarize_module.get_non_zero_cov_ROI('braf', 'V600E', cwd)
	egfr_L858R_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'L858R', cwd)
	#egfr_exon19del_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'exon19del', cwd)
	#egfr_exon20ins_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'exon20ins', cwd)
	egfr_G719X_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'G719X', cwd)
	egfr_L861Q_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'L861Q', cwd)
	egfr_S768I_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'S768I', cwd)
	egfr_T790M_cov_nonZero = summarize_module.get_non_zero_cov_ROI('egfr', 'T790M', cwd)
	#kras_G12C_cov_nonZero = summarize_module.get_non_zero_cov_ROI('kras', 'G12C', cwd)
	#kras_G13X_cov_nonZero = summarize_module.get_non_zero_cov_ROI('kras', 'G13X', cwd)
	#kras_Q61X_cov_nonZero = summarize_module.get_non_zero_cov_ROI('kras', 'Q61X', cwd)

	# fix up some of the wierd ones
	#kras_G13X_cov_nonZero['depth_gvcf'][4202] = 34
	#kras_Q61X_cov_nonZero['depth_gvcf'][6431] = 92
	#egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip('[')
	#egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip(']')
	#egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip("'")

	# fill in ROI coverage
	#summarize_module.ROI_coverage_fill_in(braf_V600E_cov_nonZero, 'BRAF', 'V600E', summaryTable)
	summarize_module.ROI_coverage_fill_in(egfr_G719X_cov_nonZero, 'EGFR', 'G719X', summaryTable)
	summarize_module.ROI_coverage_fill_in(egfr_L858R_cov_nonZero, 'EGFR', 'L858R', summaryTable)
	summarize_module.ROI_coverage_fill_in(egfr_L861Q_cov_nonZero, 'EGFR', 'L861Q', summaryTable)
	summarize_module.ROI_coverage_fill_in(egfr_S768I_cov_nonZero, 'EGFR', 'S768I', summaryTable)
	summarize_module.ROI_coverage_fill_in(egfr_T790M_cov_nonZero, 'EGFR', 'T790M', summaryTable)
	#summarize_module.ROI_coverage_fill_in(kras_G12C_cov_nonZero, 'KRAS', 'G12C', summaryTable)
	#summarize_module.ROI_coverage_fill_in(kras_G13X_cov_nonZero, 'KRAS', 'G13X', summaryTable)
	#summarize_module.ROI_coverage_fill_in(kras_Q61X_cov_nonZero, 'KRAS', 'Q61X', summaryTable)
	#summarize_module.ROI_coverage_fill_in(egfr_exon19del_cov_nonZero, 'EGFR', 'del19', summaryTable)
	#summarize_module.ROI_coverage_fill_in(egfr_exon20ins_cov_nonZero, 'EGFR', 'ins20', summaryTable)

	# trim it down
	summaryTable_trimmed = summaryTable[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'tumorCell_bool', 'mutations_found_translated']]
	summaryTable_trimmed.columns = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool', 'mutations_found']
	summaryTable_trimmed = summaryTable_trimmed[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'mutations_found', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool']]
	summaryTable_trimmed

	# add sample name col to summary table
	summaryTable_trimmed['sample_name'] = ''
	summarize_module.generic_summary_table_fill_in('sample_name', 'sample_name', summaryTable_trimmed, patientMetadata)

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
	summarize_module.validation_table_metadata_fill_in('patient_id', 'patient', validationTable_samples, patientMetadata)
	summarize_module.validation_table_metadata_fill_in('driver_gene', 'driver_gene', validationTable_samples, patientMetadata)
	summarize_module.validation_table_metadata_fill_in('driver_mutation', 'driver_mutation', validationTable_samples, patientMetadata)

	# fill in mutations found
	muts_dict = summarize_module.validation_table_dict_muts(validationTable_samples, summaryTable_trimmed)
	validationTable_samples['mutations_found'] = muts_dict.values()

	# fill in numTumorCells (various)
	tc_dict = summarize_module.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'tumorCell_bool')
	tc_cov_dict = summarize_module.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'coverage_to_ROI')
	clinMut_dict = summarize_module.validation_table_dict_generic(validationTable_samples, summaryTable_trimmed, 'clinical_mutation_found_bool')

	validationTable_samples['numTumorCells'] = tc_dict.values()
	validationTable_samples['numTumorCells_w_coverage_to_ROI'] = tc_cov_dict.values()
	validationTable_samples['numTumorCells_clinMut_found'] = clinMut_dict.values()

	# clean up 
	validationTable_samples = validationTable_samples.drop([0]) # this can change
	cols = ['sample', 'patient', 'driver_gene', 'driver_mutation', 'mutations_found', 'numTumorCells', 'numTumorCells_w_coverage_to_ROI', 'numTumorCells_clinMut_found']
	validationTable_samples = validationTable_samples[cols]

	# write
	validationTable_samples.to_csv(cwd + 'validationTable_samples.csv', index=False)

