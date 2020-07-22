''' want to test protein_variant_predictor with a larger version of
	the cosmic dataframe '''
import os
import pandas as pd
import vcfpy
import unittest
from pyfaidx import Fasta

from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.find_peptide_variants import AminoAcidMutationFinder
from cerebra.utils import *


class ProteinVariantPredictorTesterCos(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''

		data_path = os.path.abspath(__file__ + '/../' + \
										'data/test_find_peptide_variants/')
		cosmicdb_path = data_path + '/cosmic_min.tsv'
		annotation_path = data_path + '/gencode_min.gtf'
		genomefa_path = data_path + '/GRCh38_limited_chr7.fa.gz'

		self.input_path = data_path + '/vcf/'  
		self.input_paths = [self.input_path +
								x for x in os.listdir(self.input_path)]

		cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		genome_faidx = Fasta(genomefa_path)

		self.annotation_genome_tree = GenomeIntervalTree(lambda feat: feat.pos,
						(GFFFeature(row) for _,
						row in annotation_df.iterrows()))

		self.protein_variant_predictor = ProteinVariantPredictor(
									self.annotation_genome_tree, genome_faidx)

		self.aa_mutation_finder = AminoAcidMutationFinder(cosmic_df,
										annotation_df, genome_faidx, cov_bool=0)


	def test_predict_for_vcf_record(self):
		''' tests predict_for_vcf_record, get_all_overlaps and
			sequence_variants_are_equivalent methods for our test vcf set '''
		target_count = 0

		# these are all possible ensembl translation ids that should
		#	be predicted given our test vcfs
		potential_variants = ['ENSP00000415559.1:p.(Leu813Arg)',
								'ENSP00000395243.3:p.(Leu813Arg)',
								'ENSP00000275493.2:p.(Leu858Arg)',
								'ENSP00000415559.1:p.(Leu813delinsArgTrp)',
								'ENSP00000395243.3:p.(Leu813delinsArgTrp)',
								'ENSP00000275493.2:p.(Leu858delinsArgTrp)']

		for vcf_path in self.input_paths:
			curr_vcf = vcf_path.strip(self.input_path)
			vcf_reader = vcfpy.Reader.from_path(vcf_path)

			for record in vcf_reader:
				record_pos = GenomePosition.from_vcf_record(record)

				protein_variant_results = self.protein_variant_predictor \
											.predict_for_vcf_record(record)

				if not protein_variant_results:
					continue

				overlaps = self.aa_mutation_finder. \
								_cosmic_genome_tree.get_all_overlaps(record_pos)

				target_variants = (
					self.aa_mutation_finder. \
								_get_cosmic_record_protein_variant(overlap)
					for overlap in overlaps)

				for t in target_variants: # testing here buh
					print(t)

				# Filter out variants which could not be correctly obtained for
				# some reason.
				target_variants = \
						[variant for variant in target_variants if variant]

				if not target_variants:
					continue

				for result in protein_variant_results:
					predicted_variant = result.predicted_variant

					assert str(predicted_variant) in potential_variants

					if target_variants is not None:
						for target_variant in target_variants:
		
							if sequence_variants_are_equivalent(
								target_variant,
								predicted_variant,
								strict_unknown=False,
								strict_silent=True):
				
								target_count += 1
								break

		assert target_count == 1  # means that only one of these badboys 
								  #     is in COSMIC


if __name__ == "__main__":
	unittest.main()
