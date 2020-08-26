''' test cases for some of the exception handling in 
	find_peptide_variants.py '''

import os
import pandas as pd
import vcfpy
import unittest
from pathlib import Path 
from pyfaidx import Fasta

from cerebra.protein_variant_predictor import ProteinVariantPredictor
from cerebra.find_peptide_variants import AminoAcidMutationFinder
from cerebra.utils import *


class FindPeptideVariantsTester(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''

		self.data_path = os.path.abspath(__file__ + '/../' + \
									'data/test_find_peptide_variants/')

		annotation_path = self.data_path + '/gencode_min.gtf'
		genomefa_path = self.data_path + '/GRCh38_limited_chr7.fa.gz'
		cosmicdb_path = self.data_path + \
						'/CosmicGenomeScreensMutantExport.min.tsv'

		self.input_path = self.data_path + '/vcf/' 
		self.input_paths = [self.input_path +
								x for x in os.listdir(self.input_path)]

		self.cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		self.annotation_df = pd.read_csv(annotation_path, sep='\t', skiprows=5)
		self.genome_faidx = Fasta(genomefa_path)

		self.annotation_genome_tree = GenomeIntervalTree(lambda feat: 
						feat.pos,
						(GFFFeature(row) for _,
						row in self.annotation_df.iterrows()))

		self.protein_variant_predictor = ProteinVariantPredictor(
								self.annotation_genome_tree, self.genome_faidx)

		self.aa_mutation_finder = AminoAcidMutationFinder(self.cosmic_df,
									self.annotation_df, self.genome_faidx, coverage=1)


	def test_target_variant_string(self):
		''' make sure that target_variants string in 
			find_peptide_variants.py is doing what its supposed to  '''

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

				for overlap in overlaps:
					# this should be the only 'Mutation AA' str in the 
					#     COSMIC db
					assert overlap["Mutation AA"] == 'p.L858R'
					target_variant = self.aa_mutation_finder. \
							_get_cosmic_record_protein_variant(overlap)
					
					# this is the only one of our test variants
					#							 thats in COSMIC
					assert str(target_variant) == \
									'ENSP00000275493.2:p.Leu858Arg'


	def test_exception_handling(self):
		''' this case SHOULD trigger an exception in 
			find_peptide_variants.find_cell_gene_aa_mutations(). 
			explicitly feeding this function a record that we know is 
			not contained in the genome interval tree -- testing
			exception handling case '''

		for index, row in self.cosmic_df.iterrows():
			cosmic_record = self.cosmic_df.loc[index]
			mutation_aa = cosmic_record["Mutation AA"]
			gene_name = cosmic_record["Gene name"]
			transcript_accession = cosmic_record["Accession Number"]

			for tx_id, tx_record in self.protein_variant_predictor \
				.transcript_records.items():
				if tx_id.split('.')[0] == transcript_accession:
					transcript_record = tx_record
					break
			else:   # this record SHOULD trigger exception
				target_variant = self.aa_mutation_finder. \
						_get_cosmic_record_protein_variant(cosmic_record)
				
				assert target_variant == None
				assert mutation_aa != 'p.L858R'
				assert gene_name != 'EGFR' and gene_name != 'BRAF'


	def test_extract_coverage(self):
		''' testing the extract_coverage() subroutine within
		  	find_peptide_variants.find_cell_gene_aa_mutations() '''

		def init_process(curr_aa_mutation_finder):
			global current_process_aa_mutation_finder
			current_process_aa_mutation_finder = curr_aa_mutation_finder

		def process_cell(path):
			my_obj = current_process_aa_mutation_finder. \
						find_cell_gene_aa_mutations(path=path)

			return(my_obj)

		# need to rebuild this with larger version of COSMIC
		cosmicdb_path = self.data_path + '/cosmic_min.tsv'
		cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
		curr_aa_mutation_finder = AminoAcidMutationFinder(cosmic_df,
									self.annotation_df, self.genome_faidx, coverage=1)

		for vcf_path in self.input_paths:
			curr_vcf = vcf_path.strip(self.input_path)
			init_process(curr_aa_mutation_finder)
			gene_aa_mutations = process_cell(vcf_path)

			if 'A1' in curr_vcf:
				k = list(gene_aa_mutations)
				assert 'EGFR' in k

				# take a look at the variant coverage string 
				#		this is horrible -- I know
				v = gene_aa_mutations.get('EGFR')
				vl = list(v)
				ve = vl[0]
				vc = ve.split(',')[1] 
				assert vc == '[2:0]'

			elif 'A2' in curr_vcf:
				k = list(gene_aa_mutations)
				v = gene_aa_mutations.get('EGFR')
				vl = list(v)
				ve = vl[0]
				vc = ve.split(',')[1]
				assert vc == '[0:2]' or vc == '[3:2:1]'

			else:
				k = list(gene_aa_mutations)
				assert not k


if __name__ == "__main__":
	unittest.main()
