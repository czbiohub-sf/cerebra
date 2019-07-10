""" add description here """

import numpy as np
import os
from . import utils
import csv
import pandas as pd
import sys
import itertools
import warnings
import click 
from tqdm import tqdm
import multiprocessing as mp
warnings.simplefilter(action='ignore', category=FutureWarning)



""" get cmdline input """
@click.command()
@click.option('--genes_list', default = 'EGFR', prompt='path to csv file with list of genes you wish to assess coverage for', required=True, type=str)
@click.option('--nthread', default = 16, prompt='number of threads', required=True, type=int)
@click.option('--wrkdir', default = '/home/ubuntu/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def check_coverage_loci(genes_list, nthread, wrkdir):
	""" evaluate coverage to every loci in a list of user-defined genes. 
		reports on a per-loci basis. """

	print('hello world')