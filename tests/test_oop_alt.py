''' lesss try this again!! '''

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
from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.find_aa_mutations import AminoAcidMutationFinder
from cerebra.utils import *


class MyNewTestCase(unittest.TestCase):
	''' what does this test class do? '''
	@classmethod
	def setUpClass(self):
		''' set this fucker up '''
		data_path = os.path.abspath(__file__ + '/../' + 'data/test_find_aa_mutations/')
		cosmicdb_path =  data_path + '/CosmicGenomeScreensMutantExport.min.tsv'
		annotation_path = data_path + '/gencode.v33.greatestHits.annotation.gtf'
		genomefa_path = data_path + '/GRCh38.p13.genome.fa'

		input_path = data_path + '/vcf/'
		input_paths = [input_path + x for x in os.listdir(input_path)]

		self.cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		self.annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		self.genome_faidx = Fasta(genomefa_path)


	def test_aa_mutation_finder_setup(self):
		''' what are we testing here?'''
		annotation_genome_tree = GenomeIntervalTree(
			lambda feat: feat.pos,
			(GFFFeature(row) for _, row in self.annotation_df.iterrows()))

		protein_variant_predictor = ProteinVariantPredictor(
									annotation_genome_tree, self.genome_faidx)

		assert len(annotation_genome_tree.records) == 885
		assert len(annotation_genome_tree.tree_map) == 4
		assert True == True

		assert len(protein_variant_predictor.transcript_records) == 17
		assert len(protein_variant_predictor.tree.records) == 230
		assert len(protein_variant_predictor.tree.tree_map) == 2
		assert True == True


	def test_assert(self):
		assert True == True


	def test_cosmic_filter(self):
		''' what can we do here? '''
		aa_mutation_finder_bare = AminoAcidMutationFinder.__new__(AminoAcidMutationFinder)

		# this is what we want to do!
		filtered_cosmic_df = aa_mutation_finder_bare._make_filtered_cosmic_df(self.cosmic_df)

		cosmic_genome_tree = GenomeIntervalTree(
                lambda row: GenomePosition.from_str(
                    str(row["Mutation genome position"])),
                (row for _, row in filtered_cosmic_df.iterrows()))

		assert len(filtered_cosmic_df.index) == 93
		assert len(filtered_cosmic_df.columns) == 34

		assert len(cosmic_genome_tree.records) == 93
		assert len(cosmic_genome_tree.tree_map) == 1

		assert True == True


if __name__ == "__main__":
    unittest.main()

