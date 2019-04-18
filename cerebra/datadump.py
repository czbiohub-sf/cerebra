""" do s3 imports, and set up wrkdir """

import os
import click


def dir_setup():
	""" sets up the directory structure of wrkdir """
	print(' ')
	print('setting up directory structure...')

	os.system('sudo mkdir -p wrkdir')
	os.system('sudo chmod -R 777 wrkdir')

	os.system('sudo mkdir -p wrkdir/bulkVCF')
	os.system('sudo chmod -R 777 wrkdir/bulkVCF/')

	os.system('sudo mkdir -p wrkdir/coverage')
	os.system('sudo chmod -R 777 wrkdir/coverage')

	os.system('sudo mkdir -p wrkdir/filteredOut')
	os.system('sudo chmod -R 777 wrkdir/filteredOut')

	os.system('sudo mkdir -p wrkdir/fusionsOut')
	os.system('sudo chmod -R 777 wrkdir/fusionsOut')

	os.system('sudo mkdir -p wrkdir/fusions')
	os.system('sudo chmod -R 777 wrkdir/fusions')

	os.system('sudo mkdir -p wrkdir/gVCF')
	os.system('sudo chmod -R 777 wrkdir/gVCF')

	os.system('sudo mkdir -p wrkdir/scVCF')
	os.system('sudo chmod -R 777 wrkdir/scVCF')

	os.system('sudo mkdir -p wrkdir/scVCF_filtered_all')
	os.system('sudo chmod -R 777 wrkdir/scVCF_filtered_all')

	os.system('sudo mkdir -p wrkdir/test')
	os.system('sudo chmod -R 777 wrkdir/test')
	os.system('sudo mkdir -p wrkdir/test/coverage')
	os.system('sudo mkdir -p wrkdir/test/germline_filter')
	os.system('sudo mkdir -p wrkdir/test/mutationcounts_table')
	os.system('sudo mkdir -p wrkdir/test/specific_mutations')
	os.system('sudo mkdir -p wrkdir/test/summary_tables')
	os.system('sudo chmod -R 777 wrkdir/test')



""" get cmdline input """
@click.command()
@click.option('--metadata', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/metadata_all_cells_4.10.19.csv', prompt='s3 path to by-cell metadata', required=True, type=str)
@click.option('--cosmic_db', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/CosmicGenomeScreensMutantExport.tsv', prompt='s3 path to COSMIC database', required=True, type=str)
@click.option('--hg38', default = 's3://darmanis-group/singlecell_lungadeno/rawdata/hg38-plus.gtf', prompt='s3 path to hg38.gtf', required=True, type=str)
@click.option('--sc_vcf', default = 's3://lincoln.harris-work/scVCF/', prompt='s3 path to single-cell VCFs', required=True, type=str)
@click.option('--bulk_vcf', default = 's3://lincoln.harris-work/bulkVCF/', prompt='s3 path to bulk VCFs', required=True, type=str)
@click.option('--g_vcf', default = 's3://lincoln.harris-work/gVCF/', prompt='s3 path to gVCFs', required=True, type=str)



def s3_import(metadata, cosmic_db, hg38, sc_vcf, bulk_vcf, g_vcf):
	""" import necessary files from s3 """

	dir_setup()

	print('dwlding files from s3...')
	print(' ')
	
	cmd = 'aws s3 cp ' + metadata + ' wrkdir/ --quiet'
	os.system(cmd)	

	cmd = 'aws s3 cp ' + cosmic_db + ' wrkdir/ --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + hg38 + ' wrkdir/ --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + sc_vcf + ' wrkdir/scVCF/ --recursive --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + bulk_vcf + ' wrkdir/bulkVCF/ --recursive --quiet'
	os.system(cmd)

	cmd = 'aws s3 cp ' + g_vcf + ' wrkdir/gVCF/ --recursive --quiet'
	os.system(cmd)

	cmd = 'cp vcfheader.txt wrkdir/'
	os.system(cmd)

	# dwld STAR-fusion files
	cmd = 'aws s3 cp s3://lincoln.harris-work/fusions wrkdir/fusions --recursive --quiet'
	os.system(cmd)

	# dwld Seurat metadata
	cmd = 'aws s3 cp s3://lincoln.harris-work/metadataSeurat.csv wrkdir/ --quiet'
	os.system(cmd)

