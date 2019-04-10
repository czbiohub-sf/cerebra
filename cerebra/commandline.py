# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401
from functools import partial
import click
from cerebra.hello import hello
from cerebra.germlineFilter import germlineFilter

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
cli.add_command(germlineFilter, name='germlineFilter')

if __name__ == "__main__":
    cli()
