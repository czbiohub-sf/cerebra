''' really just want to run things on a big version of the cosmic df '''

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


class FindAAMutationsTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''
		data_path = os.path.abspath(__file__ + '/../' + 'data/test_find_aa_mutations/')
		cosmicdb_path =  data_path + '/cosmic_kras_egfr_braf_only.tsv.gz'   # big boi
		annotation_path = data_path + '/gencode.v33.greatestHits.annotation.gtf'
		genomefa_path = data_path + '/GRCh38_limited.fa.gz'

		self.input_path = data_path + '/vcf/'
		self.input_paths = [self.input_path + x for x in os.listdir(self.input_path)]
		
		self.cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		self.annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		self.genome_faidx = Fasta(genomefa_path)


	def test_get_all_overlaps(self):
		''' make a BARE AminoAcidMutationFinder, then setup the cosmic_genome_tree
			obj and evalutate  '''
		aa_mutation_finder_bare = AminoAcidMutationFinder.__new__(AminoAcidMutationFinder)

		# this is what we want to do!
		filtered_cosmic_df = aa_mutation_finder_bare._make_filtered_cosmic_df(self.cosmic_df)

		cosmic_genome_tree = GenomeIntervalTree(
                lambda row: GenomePosition.from_str(
                    str(row["Mutation genome position"])),
                (row for _, row in filtered_cosmic_df.iterrows()))

		# list of all possible mutations given the variants we've introduced in A1.vcf
		muts_list = ['p.L858R', 'p.G13D', 'p.G13V', 'p.G12C', 'p.G12S', 'p.G12R', 'p.G12F', 'p.G12Y']

		for vcf_path in self.input_paths:
			curr_vcf = vcf_path.strip(self.input_path)

			vcf_reader = vcfpy.Reader.from_path(vcf_path)

			for record in vcf_reader:
				record_pos = GenomePosition.from_vcf_record(record)
				overlaps = cosmic_genome_tree.get_all_overlaps(record_pos)

				for overlap in overlaps:
					mutation_aa = overlap["Mutation AA"]
		
					if curr_vcf == 'A1':
						assert str(mutation_aa) in muts_list

		assert True == True


	def test_assert(self):
		assert True == True


if __name__ == "__main__":
    unittest.main()
