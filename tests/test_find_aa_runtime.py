''' very simple tests for find_aa_mutations runtime functionality '''
import os
import pandas as pd
import filecmp
from click.testing import CliRunner
from pyfaidx import Fasta

from cerebra.find_aa_mutations import find_aa_mutations


def test_basic():
	''' does find_all_mutations returns w/o error? '''
	data_path = os.path.abspath(__file__ + '/../' +
								'data/test_find_aa_mutations/')
	vcf_path = data_path + '/vcf/A1.vcf'
	genomefa_path = data_path + '/GRCh38_limited.fa.gz'

	runner = CliRunner()
	result = runner.invoke(find_aa_mutations, [
							"--processes", 1,
							"--annotation", data_path + "/hg38-plus.min.gtf",
							"--report_coverage", 1,
							"--genomefa", genomefa_path,
							"--output", data_path + "/test_out.csv",
							vcf_path])

	assert result.exit_code == 0
	assert os.path.isfile(data_path + "/test_out.csv")

	# teardown
	os.remove(data_path + "/test_out.csv")


def test_basic_cmp():
	''' does find_all_mutations return w/o error, redux
		this one has a expected / actual file compare step '''
	from cerebra.find_aa_mutations import AminoAcidMutationFinder

	data_path = os.path.abspath(__file__ + '/../' +
									'data/test_find_aa_mutations/')
	genomefa_path = data_path + '/GRCh38_limited.fa.gz'

	annotation_path = data_path + '/gencode.v33.greatestHits.annotation.gtf'
	cov_bool = 1
	num_processes = 2
	outpath = data_path + '/test_out.csv'

	input_path = data_path + '/vcf/'
	input_paths = [input_path + x for x in os.listdir(input_path)]

	# cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
	cosmic_df = None
	annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
	genome_faidx = Fasta(genomefa_path)
	aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df,
												genome_faidx, cov_bool)

	results_df = aa_mutation_finder.find_transcript_mutations(paths=input_paths,
									processes=num_processes)

	results_df = results_df.sort_index()  # row sort
	results_df = results_df.sort_index(axis=1)  # column sort
	results_df.to_csv(outpath)

	assert os.path.isfile(outpath)

	expect_path = data_path + '/expect_out_1.csv'
	filecmp.cmp(expect_path, outpath)
