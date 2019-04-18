""" tests for cerebra::get_specific_mutations """

import os
import shutil 
from click.testing import CliRunner


def test_egfr():
	from cerebra.get_specific_mutations import get_specific_mutations

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_specific_mutations, ["--test", "True", 
													"--chrom", "7",
													"--start", "55152337", 
													"--end", "55207337", 
													"--outprefix", "egfr_test"])
	assert result.exit_code == 0
	assert os.path.isfile(wd + 'egfr_test.csv')
	assert os.path.isfile(wd + 'egfr_test_AA.csv')
	assert "egfr_L861Q,['L861Q']" in open(wd + 'egfr_test_AA.csv').read()
	assert "egfr_S768I,['S768I']" in open(wd + 'egfr_test_AA.csv').read()
	assert "egfr_L858R_test,['L858R']" in open(wd + 'egfr_test_AA.csv').read()



def test_kras():
	from cerebra.get_specific_mutations import get_specific_mutations

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_specific_mutations, ["--test", "True", 
													"--chrom", "12",
													"--start", "25209431", 
													"--end", "25250803", 
													"--outprefix", "kras_test"])
	assert result.exit_code == 0
	assert os.path.isfile(wd + 'kras_test.csv')
	assert os.path.isfile(wd + 'kras_test_AA.csv')
	assert "kras_G12C,['G12C']" in open(wd + 'kras_test_AA.csv').read()
	assert "kras_G13C,['G13C']" in open(wd + 'kras_test_AA.csv').read()



def test_braf():
	from cerebra.get_specific_mutations import get_specific_mutations

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_specific_mutations, ["--test", "True", 
													"--chrom", "7",
													"--start", "140730665", 
													"--end", "140924928", 
													"--outprefix", "braf_test"])
	assert result.exit_code == 0
	assert os.path.isfile(wd + 'braf_test.csv')
	assert os.path.isfile(wd + 'braf_test_AA.csv')
	assert "" in open(wd + 'braf_test_AA.csv').read()
	assert "braf_V600E,['V600E']" in open(wd + 'braf_test_AA.csv').read()

	