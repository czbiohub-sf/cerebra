""" do s3 imports """

import os
import click


""" get cmdline input """
@click.command()
@click.option('--metadata', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/metadata_all_cells_4.10.19.csv', prompt='s3 path to by-cell metadata', required=True, type=str)
@click.option('--cosmic_db', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/CosmicGenomeScreensMutantExport.tsv', prompt='s3 path to COSMIC database', required=True, type=str)
@click.option('--hg38', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/hg38-plus.gtf', prompt='s3 path to hg38.gtf', required=True, type=str)
@click.option('--sc_vcf', default = 's3://lincoln.harris-work/scVCF/', prompt='s3 path to single-cell VCFs', required=True, type=str)
@click.option('--bulk_vcf', default = 's3://lincoln.harris-work/bulkVCF/', prompt='s3 path to bulk VCFs', required=True, type=str)



def s3_import(metadata, cosmic_db, hg38, sc_vcf, bulk_vcf):
	""" import necessary files from s3 """
	os.system('sudo mkdir -p wrkdir')
	os.system('sudo chmod -R 777 wrkdir')
	cmd = 'aws s3 cp ' + metadata + ' wrkdir/ --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + cosmic_db + ' wrkdir/ --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + hg38 + ' wrkdir/ --quiet'
	os.system(cmd)

	os.system('sudo mkdir -p wrkdir/scVCF')
	os.system('sudo chmod -R 777 wrkdir/scVCF')
	cmd = 'aws s3 cp ' + sc_vcf + ' wrkdir/scVCF/ --recursive --quiet'
	os.system(cmd)

	os.system('sudo mkdir -p wrkdir/bulkVCF')
	os.system('sudo chmod -R 777 wrkdir/bulkVCF/')
	cmd = 'aws s3 cp ' + bulk_vcf + ' wrkdir/bulkVCF/ --recursive --quiet'
	os.system(cmd)

	cmd = 'cp vcfheader.txt wrkdir/'
	os.system(cmd)

	# only until we have coverage module written
	os.system('sudo mkdir -p wrkdir/coverage')
	os.system('sudo chmod -R 777 wrkdir/coverage')
	#cmd = 'cp /Users/lincoln.harris/code/SNP_calling_pipeline/coverage/out/*_coverageByCell.csv wrkdir/coverage'
	#os.system(cmd)

	# dwld STAR-fusion files
	os.system('sudo mkdir -p wrkdir/fusions')
	os.system('sudo chmod -R 777 wrkdir/fusions')
	cmd = 'aws s3 cp s3://lincoln.harris-work/fusions wrkdir/fusions --recursive --quiet'
	os.system(cmd)

	# dwld seurat metadata
	cmd = 'aws s3 cp s3://lincoln.harris-work/metadataSeurat.csv wrkdir/ --quiet'
	os.system(cmd)