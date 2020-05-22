  ''' basic tests for find_aa_mutations '''

import shutil 
from click.testing import CliRunner
from itertools import compress
import os
import click
import pandas as pd
import vcfpy
from pyfaidx import Fasta
import unittest

from cerebra.find_aa_mutations import find_aa_mutations
from cerebra.find_aa_mutations import AminoAcidMutationFinder
from cerebra.utils import *


class FindAAMutationsTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for FindAAMutationsTester class '''
		self.data_path = os.path.abspath(__file__ + '/../' + 'data/test_find_aa_mutations/')
		self.cosmicdb_path =  self.data_path + '/CosmicGenomeScreensMutantExport.min.tsv'
		#self.annotation_path = self.data_path + '/hg38-plus.min.gtf'
		self.annotation_path = self.data_path + '/gencode.v33.greatestHits.annotation.gtf'
	
		self.genomefa_path = self.data_path + '/GRCh38_limited.fa.gz'
		self.cov_bool = 1
		self.num_processes = 2
		self.outpath = self.data_path + '/test_out.csv'
	
		#self.input_path = '/Users/lincoln.harris/code/cerebra/cerebra/tmp/vcf/'
		self.input_path = self.data_path + '/vcf/'
		self.input_paths = [self.input_path + x for x in os.listdir(self.input_path)]

		cosmic_df = None
		annotation_df = pd.read_csv(self.annotation_path, sep='\t', skiprows=5)

		genome_faidx = Fasta(self.genomefa_path)

		self.aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df, genome_faidx, self.cov_bool)


	def test_out_basic(self):
		''' make sure we return a results_df '''
		results_df = self.aa_mutation_finder.find_transcript_mutations(paths=self.input_paths,
															processes=self.num_processes)
		results_df = results_df.sort_index() # row sort
		results_df = results_df.sort_index(axis=1) # column sort
		results_df.to_csv(self.outpath)

		assert True == True
		assert os.path.isfile(self.outpath)


	def test_out_in_depth(self):
		''' check specific attributes of results_df '''
		outfile = pd.read_csv(self.outpath, index_col=0)

		expect_index = ['A1', 'A2', 'A3', 'A4', 'A5']
		assert list(outfile.index) == expect_index

		expect_cols = ['EGFR', 'KRAS']
		assert list(outfile.columns) == expect_cols

		a1_egfr_str = outfile.loc['A1']['EGFR']
		assert 'Leu858Arg' in a1_egfr_str

		a1_kras_str = outfile.loc['A1']['KRAS']
		assert 'Gly12Ser' in a1_kras_str
		assert 'Gly13Val' in a1_kras_str

		assert True == True


	def test_helloworld(self):
		assert True == True


if __name__ == "__main__":
    unittest.main()