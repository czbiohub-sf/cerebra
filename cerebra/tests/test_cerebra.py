""" tests for cerebra """

import os
import shutil 
from click.testing import CliRunner


def test_germline_filter():
	from cerebra.germline_filter import germline_filter

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(germline_filter, ["--test", "True", "--wrkdir", wd])

	assert result.exit_code == 0
	assert os.path.isdir(wd + 'test_germline_filter/')



def test_get_specific_mutations():
	from cerebra.get_specific_mutations import get_specific_mutations

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_specific_mutations, ["--test", "True", 
													"--chrom", "7",
													"--start", "55152337", 
													"--end", "55207337", 
													"--outprefix", "testOut"])
	assert result.exit_code == 0
	assert os.path.isfile(wd + 'TEST.csv')
	assert os.path.isfile(wd + 'TEST_AA.csv')


def test_get_mutationcounts_table():
	from cerebra.get_mutationcounts_table import get_mutationcounts_table

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_mutationcounts_table, ["--nthread", "4", 
													"--test", "True", 
													"--wrkdir", wd])

	assert result.exit_code == 0
	assert os.path.isfile(wd + 'geneCellMutationCounts_artifical_TEST.csv')
	# how can i evaluate the contents of this file? 



# tear down 
wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
os.remove(wd + 'TEST.csv')
os.remove(wd + 'TEST_AA.csv')
os.remove(wd + 'testOut.csv')
os.remove(wd + 'testOut_AA.csv')
os.remove(wd + 'geneCellMutationCounts_artifical_TEST.csv')
shutil.rmtree(wd + 'test_germline_filter/')
