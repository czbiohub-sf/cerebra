# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401

# Python standard library imports
from functools import partial

# 3rd party libraries
import click

# Within-module imports
from cerebra.hello import hello


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


if __name__ == "__main__":
    cli()
