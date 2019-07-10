import io
from pathlib import Path
import random
import unittest

import pandas as pd
import numpy as np
import vcfpy
import hgvs.parser
from intervaltree import IntervalTree, Interval

from cerebra.mutation_counts import MutationCounter
from cerebra.utils import GenomePosition

class MutationCounterTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.data_path = Path(__file__).parent / "data" / "test_mutation_counts"

        self.template_vcf = (self.data_path / "template.vcf")

        self.cosmic_df = pd.read_csv((self.data_path / "CosmicGenomeScreensMutantExport.min.tsv"), sep='\t')
        self.hg38_df = pd.read_csv((self.data_path / "hg38-plus.min.gtf"), sep='\t', header=None)

        self.mutation_counter = MutationCounter(self.cosmic_df, self.hg38_df)

    def test_cosmic_subset(self):
        filtered_cosmic_df = self.mutation_counter._make_filtered_cosmic_df(self.cosmic_df)

        for _, row in filtered_cosmic_df.iterrows():
            self.assertEqual("lung", row["Primary site"])

    def test_cosmic_genome_pos_filter(self):
        lung_mut_interval_tree = IntervalTree()

        # Test for positive matches.
        for _, row in self.cosmic_df.iterrows():
            if row["Primary site"] != "lung":
                continue

            genome_pos = GenomePosition.from_str(str(row["Mutation genome position"]))

            if genome_pos is None:
                continue

            # Add the genome position to a tree for use in further assertions.
            lung_mut_interval_tree[genome_pos.start:genome_pos.end] = genome_pos.chrom

            self.assertTrue(self.mutation_counter._cosmic_subset_contains_genome_pos(genome_pos))

        # Test for negative matches, excluding negative mutation matches which
        # overlap with positive ones.
        for _, row in self.cosmic_df.iterrows():
            genome_pos = GenomePosition.from_str(str(row["Mutation genome position"]))

            if genome_pos is None:
                continue

            # genome_pos overlaps with a positive match, so it cannot be assumed
            # that it shouldn't match.
            if any(map(
                lambda it: it.data == genome_pos.chrom,
                lung_mut_interval_tree.overlap(genome_pos.start, genome_pos.end))):
                continue

            self.assertFalse(self.mutation_counter._cosmic_subset_contains_genome_pos(genome_pos))

        # Do some further negative testing to ensure that garbage genome
        # positions don't match the filter.

        negative_tests = [
            GenomePosition("nonexistent-chromosome", 0, 0),
            GenomePosition("1", -10, -1),
        ]

        for test in negative_tests:
            self.assertFalse(self.mutation_counter._cosmic_subset_contains_genome_pos(test))

    def test_gene_record_finding(self):
        for _, row in self.hg38_df.iterrows():
            genome_pos = GenomePosition.from_gtf_record(row)
            gene_record = self.mutation_counter._find_containing_gene_record(genome_pos)

            self.assertIsNotNone(gene_record)
            self.assertEqual(genome_pos, GenomePosition.from_gtf_record(gene_record))

    def test_gene_name_parsing(self):
        tests = [
            ('gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS18303";', "DDX11L1"),
            ('gene_id "MIR1302-2"; gene_name "MIR1302-2"; transcript_id "NR_036051_3"; tss_id "TSS10595";', "MIR1302-2"),
            ('gene_name "FAM138A"; transcript_id "NR_026818"; tss_id "TSS10184";', "FAM138A"),
            ('gene_name "LOC729737";', "LOC729737"),
            ('gene_name "WASH7P"', "WASH7P")
        ]

        for test, expected in tests:
            self.assertEqual(expected, self.mutation_counter._parse_gene_name(test))

    def test_mutation_count_dataframe_creation(self):
        data = [
            ("CellA", {"GENE1": 5, "GENE2": 0, "GENE3": 1            }),
            ("CellB", {"GENE1": 2,             "GENE3": 3, "GENE4": 7}),
            ("CellC", {                                              }),
        ]

        actual_df = self.mutation_counter._make_mutation_counts_df(data)
        expected_df = pd.DataFrame(index=["CellA", "CellB", "CellC"], data={
            "GENE1": np.array([5, 2, 0]),
            "GENE2": np.array([0, 0, 0]),
            "GENE3": np.array([1, 3, 0]),
            "GENE4": np.array([0, 7, 0]),
        })

        # Ensure that the indices are as expected.
        self.assertTrue(expected_df.index.equals(actual_df.index))

        # Ensure that there are no extraneous columns.
        self.assertEqual(len(expected_df), len(actual_df))

        # Ensure that each column is equal to the respective column in the
        # expected DataFrame.
        for col in list(expected_df.columns):
            self.assertTrue(expected_df[col].equals(actual_df[col]))

    @unittest.skip("""This is arguably more of an integration test than a unit
    test, and it currently fails due to inconsistencies between COSMIC and the
    reference genome.""")
    def test_cell_vcf_mutation_counting(self):
        with io.StringIO() as vcf_stream:
            with vcfpy.Reader.from_path(self.template_vcf) as template_vcf_reader:
                vcf_writer = vcfpy.Writer.from_stream(vcf_stream, header=template_vcf_reader.header)

            cosmic_subset = self.cosmic_df.loc[self.cosmic_df["Primary site"] == "lung"]

            # Write test VCF

            expected_gene_mut_counts = {}

            hgvs_parser = hgvs.parser.Parser()
            for _, row in cosmic_subset.iterrows():
                genome_pos = GenomePosition.from_str(str(row["Mutation genome position"]))

                if genome_pos is None:
                    continue

                try:
                    posedit = hgvs_parser.parse_c_posedit(row["Mutation CDS"][2:]) # pylint: disable=no-member
                except:
                    continue

                record = vcfpy.Record(
                    CHROM=genome_pos.chrom, POS=genome_pos.start + 1, ID='.', REF=posedit.edit.ref,
                    ALT=[vcfpy.Substitution(None, posedit.edit.alt)], QUAL=0, FILTER='.',
                    INFO={}
                )

                vcf_writer.write_record(record)

                gene_name = row["Gene name"]
                # Remove any gene name suffixes
                gene_name = gene_name.split('_')[0]
                expected_gene_mut_counts[gene_name] = expected_gene_mut_counts.get(gene_name, 0) + 1

            # Test mutation counting

            # Reset the buffer's cursor position
            vcf_stream.seek(0)
            _, filtered_gene_mut_counts = self.mutation_counter.find_cell_gene_mut_counts(stream=vcf_stream)

            self.assertDictEqual(expected_gene_mut_counts, filtered_gene_mut_counts)


if __name__ == "__main__":
    unittest.main()
