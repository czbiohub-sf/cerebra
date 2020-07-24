''' tests for protein_variant_predictor '''
import os
import pandas as pd
import vcfpy
import unittest
from pyfaidx import Fasta

# dont understand why this module needs to be included
from cerebra.find_peptide_variants import AminoAcidMutationFinder
from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.utils import *


class ProteinVariantPredictorTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''

		self.data_path = os.path.abspath(__file__ + '/../' + \
									'data/test_find_peptide_variants/')
		
		annotation_path = self.data_path + '/gencode_min.gtf'
		genomefa_path = self.data_path + '/GRCh38_limited_chr7.fa.gz'

		self.input_path = self.data_path + '/vcf/'
		self.input_paths = [self.input_path +
								x for x in os.listdir(self.input_path)]

		annotation_df = pd.read_csv(annotation_path, delimiter='\t', \
                                        comment='#', header=None)
		genome_faidx = Fasta(genomefa_path)

		self.annotation_genome_tree = GenomeIntervalTree(lambda feat: feat.pos,
						(GFFFeature(row) for _,
							row in annotation_df.iterrows()))

		self.protein_variant_predictor = ProteinVariantPredictor(
									self.annotation_genome_tree, genome_faidx)


	def test_protein_variant_predictor_setup(self):
		''' tests for the protein_variant_predictor class object '''
		trans_rec = self.protein_variant_predictor.transcript_records
		tree = self.protein_variant_predictor.tree.tree_map

		assert len(trans_rec) == 13
		assert len(tree) == 1

		# these are the ensembl transcript IDs that should be contained within
		#	our limited version of the cosmic db
		enst_list = ['ENST00000275493.7', 'ENST00000455089.5',
					'ENST00000342916.7', 'ENST00000454757.6',
					'ENST00000344576.6', 'ENST00000420316.6',
					'ENST00000450046.1', 'ENST00000638463.1',
					'ENST00000496384.7', 'ENST00000644969.1',
					'ENST00000646891.1', 'ENST00000288602.11',
					'ENST00000469930.2', 'ENST00000311936.8',
					'ENST00000256078.9', 'ENST00000557334.5',
					'ENST00000556131.1']

		for tx_id, tx_record in self.protein_variant_predictor \
									.transcript_records.items():
			assert tx_id in enst_list


	def test_predict_for_vcf_record(self):
		''' tests for predict_for_vcf_record, given our test vcf set '''

		# these are all possible ensembl translation ids that should be 
		#	predicted given our test vcfs

		potential_variants = ['ENSP00000415559.1:p.(Leu813Arg)',
								'ENSP00000275493.2:p.(Leu858Arg)',
								'ENSP00000395243.3:p.(Leu813Arg)',
								'ENSP00000415559.1:p.(Leu813delinsArgTrp)',
								'ENSP00000275493.2:p.(Leu858delinsArgTrp)',
								'ENSP00000395243.3:p.(Leu813delinsArgTrp)']

		for vcf_path in self.input_paths:
			vcf_reader = vcfpy.Reader.from_path(vcf_path)

			for record in vcf_reader:
				record_pos = GenomePosition.from_vcf_record(record)

				protein_variant_results = self.protein_variant_predictor \
											.predict_for_vcf_record(record)

				if not protein_variant_results:
					continue

				for result in protein_variant_results:
					predicted_variant = result.predicted_variant
					assert str(predicted_variant) in potential_variants


	def test_mito(self):
		''' TODO: description here '''

		annotation_path = self.data_path + '/gencode_min_chr7.gtf'
		genomefa_path = self.data_path + '/GRCh38_chr7_chrM.fa.gz'

		annotation_df = pd.read_csv(annotation_path, delimiter='\t', \
                                        comment='#', header=None)
		genome_faidx = Fasta(genomefa_path)

		annotation_genome_tree = GenomeIntervalTree(lambda feat: feat.pos,
						(GFFFeature(row) for _,
							row in annotation_df.iterrows()))

		protein_variant_predictor = ProteinVariantPredictor(
									annotation_genome_tree, genome_faidx)


if __name__ == "__main__":
	unittest.main()
