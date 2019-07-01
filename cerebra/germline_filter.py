"""
This program takes in two sets of vcfs, single-cell and bulk (peripheral
blood, ie. germline) and filters out the common variants found in both sc
and bulk. Creates a new directory, filteredOut, that contans the filtered vcfs.
"""

import os
import pandas as pd
import click
import vcfpy
from pathlib import Path
from tqdm import tqdm
from pathos.multiprocessing import _ProcessPool as Pool, ThreadPool
from multiprocessing import current_process, Process

from .utils import GenomePosition, GenomeIntervalTree

def create_germline_genome_tree(germline_vcf_paths):
	"""Create a `GenomeIntervalTree` comprised of the VCF records supplied by
	`germline_vcf_paths`."""

	germline_vcf_records = [list(vcfpy.Reader.from_path(path)) for path in germline_vcf_paths]

	# Flatten records
	germline_vcf_records = sum(germline_vcf_records, [])

	return GenomeIntervalTree(GenomePosition.from_vcf_record, germline_vcf_records)

def write_filtered_vcf(cell_vcf_stream, germline_tree, out_stream):
	"""Rewrite `cell_vcf_stream` to `out_stream` while removing any mutations
	found in dbSNP or shared with the VCF records in `germline_tree`."""

	cell_vcf = vcfpy.Reader.from_stream(cell_vcf_stream)
	out_vcf = vcfpy.Writer.from_stream(out_stream, header=cell_vcf.header)

	for record in cell_vcf:
		# If a record's ID field is `.`, that means that the calling software
		# did not find an ID for it in the associated database, typically dbSNP.
		# This is represented as an empty array (`[]`) in VCFPy.
		if record.ID:
			# This record is in dbSNP; skip it.
			continue

		genome_pos = GenomePosition.from_vcf_record(record)

		# Filter for only containments before checking for record equality.
		# This avoids O(n*m) complexity.
		germline_records = germline_tree.get_all_containments(genome_pos)

		for germline_record in germline_records:
			if (
				germline_record.POS == record.POS and
				germline_record.REF == record.REF and
				germline_record.ALT == record.ALT):
				break
		else:
			# Write this record if it wasn't the same as any of the containment
			# records.
			out_vcf.write_record(record)

""" launch """
@click.command()
@click.option("--processes", default=1, prompt="number of processes to use for computation", type=int)
@click.option("--germline", "germline_path", prompt="path to germline vcf files directory", required=True)
@click.option("--cells", "cells_path", prompt="path to cell vcf files directory", required=True)
@click.option("--metadata", "metadata_path", prompt="path to metadata csv file", required=True)
@click.option("--outdir", "out_path", prompt="path to output vcf files directory", required=True)
def germline_filter(processes, germline_path, cells_path, metadata_path, out_path):
	"""Given a set of single-cell VCFs and bulk VCFs (peripheral blood), this
	command removes variations from the single-cell VCFs which are shared
	with corresponding bulk VCFs or part of dbSNP."""

	germline_path = Path(germline_path)
	cells_path = Path(cells_path)
	metadata_path = Path(metadata_path)
	out_path = Path(out_path)

	metadata_df = pd.read_csv(metadata_path)

	# Create a set of all patient IDs from the metadata file.
	all_patient_ids = set(metadata_df["patient_id"])

	def process_patient(patient_id):
		# Find all non-tumor bulk VCF files for the patient ID.
		germline_wb_vcf_paths = list(germline_path.glob(patient_id + "_*_*.vcf"))

		# Fetch all cell IDs associated with the patient ID.
		cell_ids = metadata_df.loc[metadata_df["patient_id"] == patient_id]["cell_id"]

		# Use the cell IDs to create a list of all single-cell VCF files for the patient.
		cell_vcf_paths = [(cells_path / cell_id).with_suffix(".vcf") for cell_id in cell_ids]

		# Create a genome interval tree for the patient's germline bulk VCF
		# data. Only selects one germline VCF to avoid over-filtering for
		# patients with multiple germline VCFs.
		germline_tree = create_germline_genome_tree(germline_wb_vcf_paths[0:1])

		def process_cell(cell_vcf_path):
			if not cell_vcf_path.exists():
				return

			# If there were any germline VCFs for this patient, append `GF_` to
			# the file name to indicate that the output VCF was
			# germline-filtered, not just dbSNP-filtered.
			out_name_prefix = "" if len(germline_wb_vcf_paths) < 1 else "GF_"
			out_vcf_path = out_path / (out_name_prefix + cell_vcf_path.name)

			with open(cell_vcf_path, mode='r') as in_file:
				with open(out_vcf_path, mode='w') as out_file:
					write_filtered_vcf(in_file, germline_tree, out_file)

		# TODO: Remove this in Python 3.8.
		# This thread pool max-worker count is from the implementation in
		# Python 3.8. This can be removed once 3.8 is released.
		with ThreadPool(min(32, os.cpu_count() + 4)) as pool:
			pool.map(process_cell, cell_vcf_paths)

	print("Running germline filter...")
	if processes > 1:
            with Pool(processes) as pool:
                list(tqdm(pool.imap(process_patient, all_patient_ids), total=len(all_patient_ids), smoothing=0.01))
	else:
		list(map(process_patient, tqdm(all_patient_ids, smoothing=0.1)))

	print("Done!")

