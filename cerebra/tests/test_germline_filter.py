""" tests for cerebra """

from click.testing import CliRunner

def test_germline_filter():
	from cerebra.germline_filter import germline_filter

	runner = CliRunner()
	result = runner.invoke(germline_filter)

	assert result.exit_code == 0
	#assert os.path.isdir('wrkdir/test_germline_filter')
	