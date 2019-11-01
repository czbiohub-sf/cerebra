# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from cerebra.os_utils import os  # noqa: F401
from functools import partial

import click

from cerebra.germline_filter import germline_filter
from cerebra.mutation_counts import count_mutations
from cerebra.find_aa_mutations import find_aa_mutations
from cerebra.check_coverage_loci import check_coverage_loci

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])


@click.group(options_metavar='',
             subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    finds mutants in your scRNA-seq experiment
    """
    pass


cli.add_command(germline_filter, name='germline-filter')
cli.add_command(count_mutations, name='count-mutations')
cli.add_command(find_aa_mutations, name='find-aa-mutations')
cli.add_command(check_coverage_loci, name='check_coverage_loci')

if __name__ == "__main__":
    cli()
