"""
This program takes in two sets of vcfs, single-cell and bulk (peripheral
blood, ie. germline) and filters out the common variants found in both sc
and bulk. Creates a new directory, filteredOut, that contans the filtered vcfs.
"""

import pandas as pd
import click
import vcfpy
from pathlib import Path
from tqdm import tqdm
from pathos.multiprocessing import ProcessPool
from multiprocessing import current_process, Process

from .utils import GenomePosition, GenomeIntervalTree

def create_germline_genome_tree(germline_vcf_paths):
	germline_vcf_records = [list(vcfpy.Reader.from_path(path)) for path in germline_vcf_paths]

	# Flatten records
	germline_vcf_records = sum(germline_vcf_records, [])

	return GenomeIntervalTree(GenomePosition.from_vcf_record, germline_vcf_records)

def write_filtered_vcf(cell_vcf_path, germline_tree, out_stream):
	""" do the germline filter, and return a dataframe with only the
		UNIQUE entries for a given cell """

	cell_vcf = vcfpy.Reader.from_path(cell_vcf_path)
	out_vcf = vcfpy.Writer.from_stream(out_stream, header=cell_vcf.header)

	for record in cell_vcf:
		genome_pos = GenomePosition.from_vcf_record(record)
		germline_records = germline_tree.get_all_containments(genome_pos)

		for germline_record in germline_records:
			if (
				germline_record.POS == record.POS and
				germline_record.REF == record.REF and
				germline_record.ALT == record.ALT):
				break
		else:
			out_vcf.write_record(record)

""" launch """
@click.command()
@click.option("--processes", default=1, prompt="number of processes to use for computation", type=int)
@click.option("--germline", "germline_path", prompt="path to germline vcf files directory", required=True)
@click.option("--cells", "cells_path", prompt="path to cell vcf files directory", required=True)
@click.option("--metadata", "metadata_path", prompt="path to metadata csv file", required=True)
@click.option("--outdir", "out_path", prompt="path to output vcf files directory", required=True)
def germline_filter(processes, germline_path, cells_path, metadata_path, out_path):
	""" given a set of single-cell vcfs and bulk-seq vcfs (peripheral blood), this
		program subtracts out the mutations common to sc- and bulkVCF. """

	germline_path = Path(germline_path)
	cells_path = Path(cells_path)
	metadata_path = Path(metadata_path)
	out_path = Path(out_path)

	# read in patient metadata
	metadata_df = pd.read_csv(metadata_path)

	germline_vcf_paths = germline_path.glob("*.vcf")
	all_patient_ids = set((path.stem.split('_')[0] for path in germline_vcf_paths))

	def process_patient(patient_id):
		germline_wb_vcf_paths = list(germline_path.glob(patient_id + "_*_WB*.vcf"))

		cell_ids = metadata_df.loc[metadata_df["patient_id"] == patient_id]["cell_id"]
		cell_vcf_paths = [(cells_path / cell_id).with_suffix(".vcf") for cell_id in cell_ids]

		germline_tree = create_germline_genome_tree(germline_wb_vcf_paths)

		for cell_vcf_path in cell_vcf_paths:
			if not cell_vcf_path.exists():
				continue

			# If there were any whole blood VCFs for this patient, append
			# `filtered_` to the file name to indicate that the output VCF was
			# actually filtered.
			out_name_prefix = "" if len(germline_wb_vcf_paths) < 1 else "filtered_"
			out_vcf_path = out_path / (out_name_prefix + cell_vcf_path.name)

			with open(out_vcf_path, mode='w') as out_file:
				write_filtered_vcf(cell_vcf_path, germline_tree, out_file)

	print("Running germline filter...")
	if processes > 1:
            with ProcessPool(processes) as pool:
                list(tqdm(pool.imap(process_patient, all_patient_ids), total=len(all_patient_ids), smoothing=0.1))
	else:
		map(process_patient, tqdm(all_patient_ids, smoothing=0.1))

