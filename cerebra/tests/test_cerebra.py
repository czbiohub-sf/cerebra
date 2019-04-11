""" tests for cerebra """

import os
from click.testing import CliRunner

def test_germline_filter():
	from cerebra.germline_filter import germline_filter

	wd = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/'
	runner = CliRunner()
	result = runner.invoke(germline_filter, ["--test", "True", "--wrkdir", wd])

	assert result.exit_code == 0
	assert os.path.isdir(wd + 'test_germline_filter/')
	