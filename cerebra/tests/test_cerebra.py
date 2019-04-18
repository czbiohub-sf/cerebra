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
	assert os.path.isdir(wd + 'test/germline_filter/')


	
def test_get_mutationcounts_table():
	from cerebra.get_mutationcounts_table import get_mutationcounts_table

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(get_mutationcounts_table, ["--nthread", "4", 
													"--test", "True", 
													"--wrkdir", wd])

	assert result.exit_code == 0
	assert os.path.isfile(wd + 'test/mutationcounts_table/geneCellMutationCounts_artifical.csv')



def test_generate_summary_tables():
	from cerebra.generate_summary_tables import generate_summary_tables

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(generate_summary_tables, ["--test", "True", 
													"--wrkdir", wd])

	assert result.exit_code == 0
