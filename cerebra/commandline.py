# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from cerebra.os_utils import os  # noqa: F401
from functools import partial

import click

from cerebra.germline_filter import germline_filter
from cerebra.count_variants import count_variants
from cerebra.find_peptide_variants import find_peptide_variants

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])


@click.group(options_metavar='',
             subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    a tool for fast and accurate summarizing of variant calling
    format (VCF) files
    """
    pass


cli.add_command(germline_filter, name='germline-filter')
cli.add_command(count_variants, name='count-variants')
cli.add_command(find_peptide_variants, name='find-peptide-variants')

if __name__ == "__main__":
    cli()
