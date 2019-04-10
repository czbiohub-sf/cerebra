#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
# author: Lincoln
# date: 4/10/19
# script: germlineFilter.py
# 
# can we do the germline filter, for all the 'germline' variants for
# a given patient? i think we're pretty close here
# 		turning our jupyter notebook into a script
#
# trying to write this for ALL patients now
#
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
import os
import sys
import warnings

import pandas as pd
import VCF

warnings.simplefilter(action='ignore', category=FutureWarning)

#////////////////////////////////////////////////////////////////////
# writeVCF
#	wonder if we can write our own vcfs, instead of doing the csv
#   conversion
# 
#   essentially just filling in dummy vals, so that we can trick our
#     VCF reader into thinking these are properly formatted vcfs
#////////////////////////////////////////////////////////////////////
def writeVCF(df, outStr_):

	header = """##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --dbsnp /arg/3/0/hg38.vcf --native-pair-hmm-threads 16 --output /return/0 --input /arg/2/0/TH041.bam --reference /arg/0/0/hg38.fa --disable-read-filter MappingQualityReadFilter --disable-read-filter GoodCigarReadFilter --disable-read-filter NotSecondaryAlignmentReadFilter --disable-read-filter MappedReadFilter --disable-read-filter MappingQualityAvailableReadFilter --disable-read-filter NonZeroReferenceLengthAlignmentReadFilter --disable-read-filter NotDuplicateReadFilter --disable-read-filter PassesVendorQualityCheckReadFilter --disable-read-filter WellformedReadFilter  --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --max-mnp-distance 0 --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --adaptive-pruning false --do-not-recover-dangling-branches false --recover-dangling-heads false --consensus false --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --min-dangling-branch-length 4 --max-num-haplotypes-in-population 128 --min-pruning 2 --adaptive-pruning-initial-error-rate 0.001 --pruning-lod-threshold 1.0 --max-unpruned-variants 100 --debug-graph-transformations false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-use-double-precision false --debug false --use-filtered-reads-for-annotations false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --capture-assembly-failure-bam false --error-correct-reads false --do-not-run-physical-phasing false --min-base-quality-score 10 --smith-waterman JAVA --correct-overlapping-quality false --emit-ref-confidence NONE --use-new-qual-calculator true --use-old-qual-calculator false --annotate-with-num-discovered-alleles false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 30.0 --max-alternate-alleles 6 --max-genotype-count 1024 --sample-ploidy 2 --num-reference-samples-if-no-call 0 --genotyping-mode DISCOVERY --genotype-filtered-alleles false --contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false --minimum-mapping-quality 20 --disable-tool-default-annotations false --enable-all-annotations false",Version=4.1.0.0,Date="March 4, 2019 6:25:06 PM UTC">
"""
	
	df['QUAL'] = 60.28
	df['FILTER'] = '.'
	df['INFO'] = 'AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=3.00;QD=30.14;SOR=2.303'

	output_VCF = outStr_
	with open(output_VCF, 'w') as vcf:
		vcf.write(header)

	df.to_csv(output_VCF, sep="\t", mode='a', index=False)

#////////////////////////////////////////////////////////////////////
# getCellsList
#	get the list of cell names from a given patient
#////////////////////////////////////////////////////////////////////
def getPatientCellsList(scVCF_list_, patientID):
	currPatient_cells_ = []

	for item in scVCF_list_:
		currCell = item.strip('.vcf')
		currPlate = currCell.split('_')[1]
    
		rowToKeep = patientMetadata['plate'] == currPlate
    
		try:
			currPatient = patientMetadata['patient_id'][rowToKeep]
			currPatientVal = currPatient.item()

			if currPatientVal == patientID:
				currPatient_cells_.append(currCell)
		except:
			continue
	print('numCells: %d' % len(currPatient_cells_))
	return currPatient_cells_

#////////////////////////////////////////////////////////////////////
# getUniqueVCF_entries()
#     do the germline filter, and return a dataframe with only the
#     UNIQUE entries for a given cell 
#////////////////////////////////////////////////////////////////////
def getUniqueVCF_entries(patient, cell):
	basePATH = os.getcwd()
	patientPATH = basePATH + '/bulkVCF/' + patient
	cellPATH = basePATH + '/scVCF/' + cell + '.vcf'
	
	try:
		patient_df = VCF.dataframe(patientPATH)
		cell_df = VCF.dataframe(cellPATH)
	except FileNotFoundError:
		print('FILE NOT FOUND: %s' % cellPATH)
		return
    
	patient_df_trimmed = patient_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
	cell_df_trimmed = cell_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
    
	# get whats SHARED between patient and cell 
	#    FIND GERMLINE MUTATIONS
	patient_cell_concat = pd.concat([patient_df_trimmed, cell_df_trimmed])
	rowsToKeep = patient_cell_concat.duplicated()
	patient_cell_shared = patient_cell_concat[rowsToKeep]
	patient_cell_shared = patient_cell_shared.reset_index(drop=True)

	# now go back to the original cell df, pull out whats UNIQUE 
	#     THIS IS THE GERMLINE FILTER!!
	cell_cell_concat = pd.concat([cell_df_trimmed, patient_cell_shared])
	cell_cell_concat_noDups = cell_cell_concat.drop_duplicates(keep=False)
	cell_cell_concat_noDups = cell_cell_concat_noDups.reset_index(drop=True)
    
	return(cell_cell_concat_noDups)

#////////////////////////////////////////////////////////////////////
# main()
#	get the patient name from the cmdline, set up patientMetadata, 
#   call subroutines to get list of per-patient cells, then call the 
#   filtering func and write output to new csv
#////////////////////////////////////////////////////////////////////

global patientMetadata

# read in patient metadata
patientMetadata = pd.read_csv('../metadata_all_cells_4.10.19.csv')

# get a list of all the single-cell VCF files
cwd = os.getcwd()
vcfDir = cwd + '/scVCF/'
scVCF_list = os.listdir(vcfDir)

# get list of bulk VCF files
bulkVCF_dir = cwd + '/bulkVCF/'
bulkVCF_list = os.listdir(bulkVCF_dir)

patientsRun = [] # need to keep track of which patients have been run

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
		currPatient_cells = getPatientCellsList(scVCF_list, currPatient)

		# inner loop -- by CELL 
		for currCell in currPatient_cells:
			currCell_unique = getUniqueVCF_entries(item, currCell)
			outStr = './filteredOut/' + currCell + '_unique.vcf'
			#currCell_unique.to_csv(outStr, index=False)
			writeVCF(currCell_unique, outStr)
			#continue
			
		patientsRun.append(currPatient)

#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////

@click.command()
@click.option('--count', default=5, help='Number of greetings.')
@click.option('--name', prompt='Your name',
              help='The person to greet.')

def germlineFilter(count, name):
    """Simple program that greets NAME for a total of COUNT times, in color."""
    for x in tqdm(range(count)):
        # note that colorama.init() doesn't need to be called for the colors
        # to work
        click.echo(click.style('Hello %s!' % name, fg=random.choice(COLORS)))
