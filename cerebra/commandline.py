# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401
from functools import partial

import click

from cerebra.hello import hello
from cerebra.germline_filter import germline_filter
from cerebra.get_mutationcounts_table import get_mutationcounts_table

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])

@click.group(options_metavar='', subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    finds mutants in your scRNA-seq experiment

    not sure what this func does. must be an __init__
    sort of thing? 
    """
    pass

cli.add_command(hello, name='hello')
cli.add_command(germline_filter, name='germline_filter')
cli.add_command(get_mutationcounts_table, name='get_mutationcounts_table')

if __name__ == "__main__":
    cli()
