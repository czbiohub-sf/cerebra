import unittest
import math
import io
from pathlib import Path

import vcfpy

from cerebra.germline_filter import write_filtered_vcf
from cerebra.utils import GenomePosition, GenomeIntervalTree

class GermlineFilterTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.data_path = Path(__file__).parent / "data" / "test_germline_filter"

    def test_germline_filter(self):


        filtered_cell_vcf_paths = self.data_path.glob("GF_*.vcf")

        for filtered_cell_vcf_path in filtered_cell_vcf_paths:
            cell_name = filtered_cell_vcf_path.stem.replace("GF_", '')
            cell_vcf_path = self.data_path / (cell_name + ".vcf")
            germline_vcf_paths = self.data_path.glob(cell_name + "_GL*.vcf")

            # Create germline genome tree

            germline_vcf_records = [list(vcfpy.Reader.from_path(path)) for path in germline_vcf_paths]

            # Flatten records
            germline_vcf_records = sum(germline_vcf_records, [])

            germline_genome_tree = GenomeIntervalTree(GenomePosition.from_vcf_record, germline_vcf_records)

            # Test writing VCF

            with io.StringIO() as out_file:
                with open(cell_vcf_path, mode='r') as in_file:
                    write_filtered_vcf(in_file, germline_genome_tree, out_file)

                # Reset the buffer's cursor position
                out_file.seek(0)

                with open(filtered_cell_vcf_path, mode='r') as expected_file:
                    expected_reader = vcfpy.Reader.from_stream(expected_file)
                    out_reader = vcfpy.Reader.from_stream(out_file)

                    expected_records, out_records = list(expected_reader), list(out_reader)

                    self.assertEqual(expected_records, out_records)


if __name__ == "__main__":
    unittest.main()
