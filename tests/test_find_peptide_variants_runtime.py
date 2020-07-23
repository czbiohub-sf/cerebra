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
	genomefa_path = data_path + '/GRCh38_limited_chr7.fa.gz'
	cosmicdb_path = data_path + '/CosmicGenomeScreensMutantExport.min.tsv'

	runner = CliRunner()
	result = runner.invoke(find_peptide_variants, [
							"--processes", 1,
							"--annotation", data_path + "/gencode_min.gtf",
							"--cosmicdb", cosmicdb_path,
							"--report_coverage", 0,
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
	genomefa_path = data_path + '/GRCh38_limited_chr7.fa.gz'

	annotation_path = data_path + '/gencode_min.gtf'
	cov_bool = 0
	num_processes = 2   # want to include that multiprocessing module
	outpath = data_path + '/test_out.csv'

	input_path = data_path + '/vcf/'
	input_paths = [input_path + x for x in os.listdir(input_path)]

	cosmic_df = None
	annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
	genome_faidx = Fasta(genomefa_path)
	aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df,
												genome_faidx, cov_bool)

	results_df = aa_mutation_finder.find_transcript_mutations(
									paths=input_paths,
									processes=num_processes)

	results_df = results_df.sort_index()  # row sort
	results_df = results_df.sort_index(axis=1)  # column sort
	results_df.to_csv(outpath)

	assert os.path.isfile(outpath)

	expect_path = data_path + '/expect_out.csv'
	filecmp.cmp(expect_path, outpath)

	# teardown
	os.remove(data_path + "/test_out.csv")
	#os.remove(data_path + "/test_out.json")
