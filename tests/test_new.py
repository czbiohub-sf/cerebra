''' lets git it buh '''

import os
import pandas as pd
import vcfpy
import unittest
from pyfaidx import Fasta

from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.find_peptide_variants import AminoAcidMutationFinder
from cerebra.utils import *

class PeptideVariantsTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''

		data_path = os.path.abspath(__file__ + '/../' + \
									'data/test_find_peptide_variants/')

		annotation_path = data_path + '/gencode_min.gtf'
		genomefa_path = data_path + '/GRCh38_limited_chr7.fa.gz'
		cosmicdb_path = data_path + '/cosmic_min.tsv'

		self.input_path = data_path + '/vcf/' 
		self.input_paths = [self.input_path +
								x for x in os.listdir(self.input_path)]

		cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		genome_faidx = Fasta(genomefa_path)

		self.annotation_genome_tree = GenomeIntervalTree(lambda feat: 
						feat.pos,
						(GFFFeature(row) for _,
						row in annotation_df.iterrows()))

		self.protein_variant_predictor = ProteinVariantPredictor(
								self.annotation_genome_tree, genome_faidx)

		self.aa_mutation_finder = AminoAcidMutationFinder(cosmic_df,
									annotation_df, genome_faidx, cov_bool=0)


	def test_helloworld(self):
		assert True == True


	def test_scratchpad(self):
		''' what we doing here buh? '''

		out = self.aa_mutation_finder.find_transcript_mutations(
										self.input_paths, processes=2)

		assert True == True


if __name__ == "__main__":
	unittest.main()
