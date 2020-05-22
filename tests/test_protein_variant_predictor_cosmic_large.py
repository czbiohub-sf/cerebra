''' want to increase test cov in protein_variant_predictor '''

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


class ProteinVariantPredictorTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''

		data_path = os.path.abspath(__file__ + '/../' + 'data/test_find_aa_mutations/')
		cosmicdb_path =  data_path + '/cosmic_kras_egfr_braf_only.tsv.gz'
		annotation_path = data_path + '/gencode.v33.greatestHits.annotation.gtf'
		genomefa_path = data_path + '/GRCh38_limited.fa.gz'

		self.input_path = data_path + '/vcf/'
		self.input_paths = [self.input_path + x for x in os.listdir(self.input_path)]
		
		cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		genome_faidx = Fasta(genomefa_path)

		self.annotation_genome_tree = GenomeIntervalTree(lambda feat: feat.pos, 
						(GFFFeature(row) for _, row in annotation_df.iterrows()))

		self.protein_variant_predictor = ProteinVariantPredictor(
									self.annotation_genome_tree, genome_faidx)

		self.aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df, 
															genome_faidx, cov_bool=0)


	def test_predict_for_vcf_record(self):
		''' todo: add description '''
		target_count = 0
		potential_variants = ['ENSP00000415559.1:p.(Leu813Arg)', 
								'ENSP00000395243.3:p.(Leu813Arg)', 
								'ENSP00000275493.2:p.(Leu858Arg)', 
								'ENSP00000308495.3:p.(Gly12Ser)', 
								'ENSP00000452512.1:p.(Gly12Ser)', 
								'ENSP00000451856.1:p.(Gly12Ser)', 
								'ENSP00000256078.4:p.(Gly12Ser)', 
								'ENSP00000308495.3:p.(Gly13Val)', 
								'ENSP00000452512.1:p.(Gly13Val)', 
								'ENSP00000451856.1:p.(Gly13Val)', 
								'ENSP00000256078.4:p.(Gly13Val)', 
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

				overlaps = self.aa_mutation_finder._cosmic_genome_tree.get_all_overlaps(record_pos)

				target_variants = (
					self.aa_mutation_finder._get_cosmic_record_protein_variant(overlap)
					for overlap in overlaps)

				# Filter out variants which could not be correctly obtained for
				# some reason.
				target_variants = [variant for variant in target_variants if variant]

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

		assert target_count == 3
		assert True == True


	def test_assert(self):
		assert True == True


if __name__ == "__main__":
    unittest.main()