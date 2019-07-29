""" tests for check_coverage_loci """

import os
import filecmp
from click.testing import CliRunner


def test_cov_module():
	from cerebra.check_coverage_loci import check_coverage_loci

	genes_file = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/genesList.csv'
	n_thread = 2
	prefix = 'test_cov_mod_out.csv'
	vcf_dir_ = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/vcf_test_set/'
	cosmic_db_ = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/CosmicGenomeScreensMutantExport.tsv'

	runner = CliRunner()
	result = runner.invoke(check_coverage_loci, ["--genes_list", genes_file,
						 "--nthread", n_thread, "--outprefix", prefix, 
						 "--vcf_dir", vcf_dir_, "--cosmic_db", cosmic_db_])

	assert result.exit_code == 0
	assert os.path.isfile(prefix)

	p2 = 'artificial_vcf_gold_std.csv'
	assert filecmp.cmp(prefix,p2), 'filecmp to artificial_vcf_gold_std.csv FAILED'