""" tests for check_coverage_loci """

import os
import filecmp
from click.testing import CliRunner


def test_cov_module():
	from cerebra.check_coverage_loci import check_coverage_loci

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	genes_file = 'genesList.csv'
	n_thread = 2
	prefix = 'test_cov_mod_out.csv'

	runner = CliRunner()
	result = runner.invoke(check_coverage_loci, ["--genes_list", genes_file,
						 "--nthread", n_thread, "--outprefix", prefix, 
						 "--wrkdir", wd])

	assert result.exit_code == 0
	assert os.path.isfile(wd + prefix)

	p1 = wd + prefix
	p2 = '/Users/lincoln.harris/Desktop/artificial_vcf_gold_std.csv'
	assert filecmp.cmp(p1,p2)