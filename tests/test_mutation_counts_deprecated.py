""" tests for mutation_counts module """

import os
import filecmp
from click.testing import CliRunner


def test_mutations_counts():
	from cerebra.mutation_counts import mutation_counts

	n_thread = 16
	prefix = 'test_mutation_counts_out.csv'
	vcf_dir_ = '/home/lincoln/code/cerebra/cerebra/wrkdir/vcf_test_set/*'
	cosmic_db_ = '/home/lincoln/code/cerebra/cerebra/wrkdir/CosmicGenomeScreensMutantExport.csv'
	hg38_anno = '/home/lincoln/code/cerebra/cerebra/wrkdir/hg38-plus.gtf'

	runner = CliRunner()
	result = runner.invoke(mutation_counts, ["--processes", n_thread,
						 "--cosmicdb", cosmic_db_, "--refgenome", hg38_anno, 
						 "--outfile", prefix, vcf_dir_])

	assert result.exit_code == 0
	assert os.path.isfile(prefix)

	p2 = 'artificial_vcf_gold_std_counts.csv'
	assert filecmp.cmp(prefix,p2), 'filecmp to artificial_vcf_gold_std_counts.csv FAILED'