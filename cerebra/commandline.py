# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401
from functools import partial

import click

from cerebra.hello import hello
from cerebra.datadump import s3_import
from cerebra.germline_filter import germline_filter
from cerebra.get_mutationcounts_table import get_mutationcounts_table
from cerebra.get_aa_mutations import get_aa_mutations
from cerebra.get_mutationalburden import get_mutationalburden
from cerebra.generate_summary_tables import generate_summary_tables
from cerebra.generate_summary_tables_test import generate_summary_tables_test
from cerebra.fusion_search import fusion_search
from cerebra.fusions_x_cell import fusions_x_cell
from cerebra.check_coverage import check_coverage

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])

@click.group(options_metavar='', subcommand_metavar='<command>',
             context_settings=settings)

def cli():
    """
    finds mutants in your scRNA-seq experiment
    """
    pass

cli.add_command(hello, name='hello')
cli.add_command(s3_import, name='s3_import')
cli.add_command(germline_filter, name='germline_filter')
cli.add_command(get_mutationcounts_table, name='get_mutationcounts_table')
cli.add_command(get_aa_mutations, name='get_aa_mutations')
cli.add_command(get_mutationalburden, name='get_mutationalburden')
cli.add_command(generate_summary_tables, name='generate_summary_tables')
cli.add_command(generate_summary_tables_test, name='generate_summary_tables_test')
cli.add_command(fusion_search, name='fusion_search')
cli.add_command(fusions_x_cell, name='fusions_x_cell')
cli.add_command(check_coverage, name='check_coverage')

if __name__ == "__main__":
    cli()
