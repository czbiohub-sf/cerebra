''' adnl tests for mutation_counts '''
import unittest
import os
import pandas as pd
import vcfpy
from click.testing import CliRunner

from cerebra.count_variants import MutationCounter
from cerebra.utils import *


class TestMutationCounter(unittest.TestCase):
	@classmethod
	def setUpClass(self):
		''' __init__ method for class obj '''
		self.data_path = os.path.abspath(__file__ + '/../' + \
											'data/test_variant_counts/')
		self.cosmicdb_path = self.data_path + \
										'/cosmic_kras_egfr_braf_only.tsv.gz'
		self.annotation_path = self.data_path + '/hg38-plus.sub.gtf'

		self.input_path = self.data_path + '/vcf/'
		self.input_paths = [self.input_path + \
								x for x in os.listdir(self.input_path)]

		self.cosmic_df = pd.read_csv(self.cosmicdb_path, delimiter='\t')
		self.refgenome_df = pd.read_csv(self.annotation_path,
											delimiter='\t', header=None)


	def test_mutation_counter_init(self):
		''' make a bare MutationCounter object and make sure basic attributes
			set up correctly '''

		mutation_counter = MutationCounter.__new__(MutationCounter)
		filtered_cosmic_df = mutation_counter. \
								_make_filtered_cosmic_df(self.cosmic_df)
		cosmic_genome_tree = GenomeIntervalTree(
								lambda row: GenomePosition.from_str(
								str(row["Mutation genome position"])),
								(record for idx, record in
									filtered_cosmic_df.iterrows()))

		hg38_genome_tree = GenomeIntervalTree(
								GenomePosition.from_gtf_record,
								(record for idx, record in
									self.refgenome_df.iterrows()))

		assert len(cosmic_genome_tree.tree_map) == 2
		assert len(cosmic_genome_tree.records) == 544

		assert len(hg38_genome_tree.tree_map) == 2
		assert len(hg38_genome_tree.records) == 211


	def test_find_cell_gene_mut_counts(self):
		''' set up a bare MutationCounter object with all of the correct
			attributes, then test to see if find_cell_gene_mut_counts returns
			correctly, for our test vcf set '''

		mutation_counter = MutationCounter.__new__(MutationCounter)

		filtered_cosmic_df = mutation_counter. \
								_make_filtered_cosmic_df(self.cosmic_df)
		mutation_counter._cosmic_genome_tree = GenomeIntervalTree(
								lambda row: GenomePosition.from_str(
								str(row["Mutation genome position"])),
								(record for idx,
									record in filtered_cosmic_df.iterrows()))
  
		mutation_counter._hg38_genome_tree = GenomeIntervalTree(
									GenomePosition.from_gtf_record,
									(record for idx,
										record in self.refgenome_df.iterrows()))

		a1_expect = {'EGFR': 2, 'KRAS': 2}

		for vcf in self.input_paths:
			curr_vcf = vcf.strip(self.input_path)
			counts = mutation_counter.find_cell_gene_mut_counts(path=vcf)

			if 'A1' in curr_vcf:
				assert counts[0] == a1_expect
			else:
				assert counts[0] == {}


	def test_parse(self):
		''' set up a bare MutationCounter object with all of the correct
			attributes, then test to see if _parse_gene_names returns
			correctly, for our test vcf set '''

		mutation_counter = MutationCounter.__new__(MutationCounter)

		mutation_counter._hg38_genome_tree = GenomeIntervalTree(
									GenomePosition.from_gtf_record,
									(record for idx,
										record in self.refgenome_df.iterrows()))

		A1_gene_names = ['KRAS', 'EGFR']

		for vcf in self.input_paths:
			curr_vcf = vcf.strip(self.input_path)
			vcf_reader = vcfpy.Reader.from_path(vcf)

			for record in vcf_reader:
				genome_pos = GenomePosition.from_vcf_record(record)
				gene_row = mutation_counter. \
								_find_containing_gene_record(genome_pos)

				if gene_row is None:
					continue

				gene_name = mutation_counter._parse_gene_name(gene_row[8])

				if 'A1' in curr_vcf:
					assert gene_name in A1_gene_names
				else:
					assert gene_name == ''


	def test_runtime(self):
		''' full runtime test
			does count_variants returns w/o error? '''
		from cerebra.count_variants import count_variants

		runner = CliRunner()
		result = runner.invoke(count_variants, ["--processes", 1,
										"--cosmicdb", self.cosmicdb_path,
										"--refgenome", self.annotation_path,
										"--outfile", self.data_path +
										'/test_out.csv',
										self.input_path + 'A1.vcf'])

		assert result.exit_code == 0
		assert os.path.isfile(self.data_path + "/test_out.csv")

		# teardown
		os.remove(self.data_path + "/test_out.csv")


if __name__ == "__main__":
	unittest.main()
