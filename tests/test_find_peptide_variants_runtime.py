''' very simple tests for find_peptide_variants runtime functionality '''
import os
import pandas as pd
import filecmp
from click.testing import CliRunner
from pyfaidx import Fasta

from cerebra.find_peptide_variants import find_peptide_variants


def test_basic():
	''' does find_all_mutations returns w/o error? 
		basic runtime test '''
	data_path = os.path.abspath(__file__ + '/../' +
								'data/test_find_peptide_variants/')
	vcf_path = data_path + '/vcf/A1.vcf'
	cosmicdb_path = data_path + '/CosmicGenomeScreensMutantExport.min.tsv'

	genomefa_gz_path = data_path + '/GRCh38_limited_chr7.fa.gz'
	genomefa_path = data_path + '/GRCh38_limited_chr7.fa'
	cmd = 'gunzip -c ' + genomefa_gz_path + ' > ' + genomefa_path
	os.system(cmd)

	runner = CliRunner()
	result = runner.invoke(find_peptide_variants, [
							"--processes", 2,
							"--annotation", data_path + "/gencode_min.gtf",
							"--cosmicdb", cosmicdb_path,
							"--genomefa", genomefa_path,
							"--output_path", data_path + "/test_out.csv",
							vcf_path])

	assert result.exit_code == 0
	assert os.path.isfile(data_path + "/test_out.csv")
	assert os.path.isfile(data_path + "/test_out.json")

	# teardown
	os.remove(data_path + "/test_out.csv")
	os.remove(data_path + "/test_out.json")


def test_gtf():
	''' lets try the same thing with a different gtf '''
	data_path = os.path.abspath(__file__ + '/../' +
								'data/test_find_peptide_variants/')
	vcf_path = data_path + '/vcf/A1.vcf'
	genomefa_path = data_path + '/GRCh38_limited_chr7.fa'
	cosmicdb_path = None

	runner = CliRunner()
	result = runner.invoke(find_peptide_variants, [
							"--processes", 2,
							"--annotation", data_path + "/hg38-plus.min.gtf",
							"--cosmicdb", cosmicdb_path,
							"--genomefa", genomefa_path,
							"--output_path", data_path + "/test_out.csv",
							vcf_path])

	assert result.exit_code == 0
	assert os.path.isfile(data_path + "/test_out.csv")
	assert os.path.isfile(data_path + "/test_out.json")

	# teardown
	os.remove(data_path + "/test_out.csv")
	os.remove(data_path + "/test_out.json")


def test_basic_cmp():
	''' does find_all_mutations return w/o error, redux
		this one has an expected vs. actual file compare step 
		pytest keeps telling me filecmp.cmp() is outdated? '''
	from cerebra.find_peptide_variants import AminoAcidMutationFinder

	data_path = os.path.abspath(__file__ + '/../' +
									'data/test_find_peptide_variants/')
	genomefa_path = data_path + '/GRCh38_limited_chr7.fa'

	annotation_path = data_path + '/gencode_min.gtf'
	coverage = 0
	num_processes = 1   # want to include that multiprocessing module
	outpath = data_path + '/test_out.csv'

	input_path = data_path + '/vcf/'
	input_paths = [input_path + x for x in os.listdir(input_path)]

	cosmic_df = None
	annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
	genome_faidx = Fasta(genomefa_path)
	aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df,
												genome_faidx, coverage)

	results_df = aa_mutation_finder.find_transcript_mutations(
									paths=input_paths,
									processes=num_processes)

	results_df = results_df.sort_index()  # row sort
	results_df = results_df.sort_index(axis=1)  # column sort
	results_df.to_csv(outpath)

	assert os.path.isfile(outpath)

	outdf = pd.read_csv(outpath, index_col=0)
	expect_index = ['A1', 'A2', 'A3', 'A4', 'A5']
	expect_cols = ['EGFR']

	assert list(outdf.index) == expect_index
	assert list(outdf.columns) == expect_cols

	indel_ensps = ['ENSP00000275493.2:p.(Leu858delinsArgTrp)', 
					'ENSP00000395243.3:p.(Leu813delinsArgTrp)', 
					'ENSP00000415559.1:p.(Leu813delinsArgTrp)']

	snp_ensps = ['ENSP00000415559.1:p.(Leu813Arg)', 
					'ENSP00000395243.3:p.(Leu813Arg)', 
					'ENSP00000275493.2:p.(Leu858Arg)', ]

	a1_matches = outdf.loc['A1']['EGFR']
	a2_matches = outdf.loc['A2']['EGFR']

	# this is a bit backwards but ok considering the wierd format of outfile
	for x in indel_ensps:
		assert x in a1_matches

	for y in snp_ensps:
		assert y in a1_matches

	for z in snp_ensps:
		assert z in a2_matches

	#teardown
	os.remove(data_path + "/test_out.csv")
	#os.remove(data_path + "/test_out.json")
