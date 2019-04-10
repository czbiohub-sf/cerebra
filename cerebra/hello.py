"""
This is an example subcommand for a git-like interface

Cribbed from the Click documentation https://click.palletsprojects.com/en/7.x/
"""
import random

import click
from tqdm import tqdm

# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.
from .os_utils import os    # noqa: F401


COLORS = 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan'


@click.command()
@click.option('--count', default=5, help='Number of greetings.')
@click.option('--name', prompt='Your name',
              help='The person to greet.')

def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times, in color."""
    for x in tqdm(range(count)):
        # note that colorama.init() doesn't need to be called for the colors
        # to work
        click.echo(click.style('Hello %s!' % name, fg=random.choice(COLORS)))
