''' testing some of the key modules of find_peptide_variants '''
import os
import pandas as pd
import vcfpy
import unittest
from pyfaidx import Fasta

from cerebra.find_peptide_variants import AminoAcidMutationFinder
from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.utils import *


class FindAAMutationsTesterMod(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''
		data_path = os.path.abspath(__file__ + '/../' +
										'data/test_find_peptide_variants/')
		cosmicdb_path = data_path + '/cosmic_min.tsv'

		annotation_path = data_path + '/gencode_min.gtf'
		genomefa_path = data_path + '/GRCh38_limited_chr7.fa.gz'

		self.input_path = data_path + '/vcf/'
		self.input_paths = [self.input_path + x for x in os.listdir(self.input_path)]

		self.cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		self.annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		self.genome_faidx = Fasta(genomefa_path)


	def test_aa_mutation_finder_setup(self):
		''' essentially making an AminoAcidMutationFinder object, from scratch '''
		annotation_genome_tree = GenomeIntervalTree(
			lambda feat: feat.pos,
			(GFFFeature(row) for _, row in self.annotation_df.iterrows()))

		protein_variant_predictor = ProteinVariantPredictor(
									annotation_genome_tree, self.genome_faidx)

		assert len(annotation_genome_tree.records) == 827
		assert len(annotation_genome_tree.tree_map) == 2

		assert len(protein_variant_predictor.transcript_records) == 13
		assert len(protein_variant_predictor.tree.records) == 218
		assert len(protein_variant_predictor.tree.tree_map) == 1


	def test_cosmic_filter(self):
		''' make a BARE AminoAcidMutationFinder, then setup the cosmic_genome_tree
			obj and evalutate  '''
		aa_mutation_finder_bare = AminoAcidMutationFinder. \
									__new__(AminoAcidMutationFinder)

		filtered_cosmic_df = aa_mutation_finder_bare. \
									_make_filtered_cosmic_df(self.cosmic_df)

		cosmic_genome_tree = GenomeIntervalTree(
								lambda row: GenomePosition.
								from_str(str(row["Mutation genome position"])),
								(row for _, row in filtered_cosmic_df.iterrows()))


		assert len(filtered_cosmic_df.index) == 283
		assert len(filtered_cosmic_df.columns) == 34

		assert len(cosmic_genome_tree.records) == 283
		assert len(cosmic_genome_tree.tree_map) == 1


	def test_genome_position_init(self):
		''' testing genome pos strings to make sure they contain what we've
			placed in the test vcf files '''

		A1_gps = ['7:55191822-55191822', '7:55191822-55191822']
		A2_gps = ['1:631862-631862', '1:633561-633561', '1:634112-634112',
					 '1:634229-634229', '1:914949-914949',
					 '7:55191822-55191822']
		A3_gps = ['1:629906-629906', '1:634112-634112', '1:634229-634229',
					'1:634244-634244', '1:1013541-1013541']
		A4_gps = ['1:1010878-1010878', '1:1010892-1010892',
					'1:1010895-1010895', '1:1010904-1010904',
					'1:1010905-1010905']
		A5_gps = ['1:630317-630317', '1:633747-633747', '1:633824-633824',
					'1:1233917-1233917', '1:1816741-1816741']

		for vcf_path in self.input_paths:
			curr_vcf = vcf_path.strip(self.input_path)

			vcf_reader = vcfpy.Reader.from_path(vcf_path)
			curr_gps = []

			for record in vcf_reader:
				record_pos = GenomePosition.from_vcf_record(record)
				curr_gps.append(str(record_pos))

			if curr_vcf == 'A1':
				assert set(curr_gps) == set(A1_gps)
			elif curr_vcf == 'A2':
				assert set(curr_gps) == set(A2_gps)
			elif curr_vcf == 'A3':
				assert set(curr_gps) == set(A3_gps)
			elif curr_vcf == 'A4':
				assert set(curr_gps) == set(A4_gps)
			elif curr_vcf == 'A5':
				assert set(curr_gps) == set(A5_gps)


if __name__ == "__main__":
	unittest.main()
