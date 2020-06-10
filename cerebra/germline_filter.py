import os
import pandas as pd
import click
import vcfpy
from pathlib import Path
from tqdm import tqdm
from pathos.multiprocessing import _ProcessPool as Pool, ThreadPool

from .utils import GenomePosition, GenomeIntervalTree


def create_germline_genome_tree(germline_vcf_paths):
    """Create a `GenomeIntervalTree` comprised of the VCF records supplied by
    `germline_vcf_paths`."""

    germline_vcf_records = [list(vcfpy.Reader.from_path(path)) for path in
                            germline_vcf_paths]

    # Flatten records
    germline_vcf_records = sum(germline_vcf_records, [])

    return GenomeIntervalTree(GenomePosition.from_vcf_record,
                              germline_vcf_records)


def write_filtered_vcf(cell_vcf_stream, germline_tree, out_stream):
    """Rewrite `cell_vcf_stream` to `out_stream` while removing any mutations
    found in dbSNP or shared with the VCF records in `germline_tree`."""

    cell_vcf = vcfpy.Reader.from_stream(cell_vcf_stream)
    out_vcf = vcfpy.Writer.from_stream(out_stream, header=cell_vcf.header)

    for record in cell_vcf:
        # If a record's ID field is `.`, that means that the calling
        # software did not find an ID for it in the associated database,
        # typically dbSNP. This is represented as an empty array (`[]`)
        # in VCFPy.

        if record.ID:   # This record is in dbSNP; skip it.
            continue

        genome_pos = GenomePosition.from_vcf_record(record)

        # Filter for only containments before checking for record equality.
        # This avoids O(n*m) complexity.
        germline_records = germline_tree.get_all_containments(genome_pos)

        unique_alts = set(record.ALT)

        for germline_record in germline_records:
            if (
                germline_record.POS != record.POS or
                germline_record.REF != record.REF):
                continue

            unique_alts -= set(germline_record.ALT)

        # Write this record if there are still has remaining (unique) ALTs.
        if unique_alts:
            record.ALT = list(unique_alts)
            out_vcf.write_record(record)


""" launch """


@click.command()
@click.option("--processes", default=1,
              prompt="number of processes to use for computation", type=int)
@click.option("--control_path", "control_path",
              prompt="path to control/germline samples vcf files directory", required=True)
@click.option("--experimental_path", "experimental_path",
              prompt="path to experimental samples vcf files directory", required=True)
@click.option("--metadata", "metadata_path",
              prompt="path to metadata csv file", required=True)
@click.option("--outdir", "out_path",
              prompt="path to output vcf files directory", required=True)
def germline_filter(processes, control_path, experimental_path, metadata_path,
                    out_path):
    """ filter out common SNPs/indels between control/germline samples and samples
        of interest """
    control_path = Path(control_path)
    experimental_path = Path(experimental_path)
    metadata_path = Path(metadata_path)
    out_path = Path(out_path)

    metadata_df = pd.read_csv(metadata_path)

    # Create a set of all patient IDs from the metadata file.
    all_germline_sample_ids = set(metadata_df["germline_sample_id"])

    def process_patient(germline_sample_id):
        # Find all non-tumor bulk VCF files for the patient ID.
        germline_wb_vcf_paths = list(
            control_path.glob(germline_sample_id + "_*_*.vcf"))

        # Fetch all cell IDs associated with the patient ID.
        experimental_sample_ids = metadata_df.loc[
                                metadata_df["germline_sample_id"] == 
                                germline_sample_id]["experimental_sample_id"]

        # Use the cell IDs to create a list of all single-cell VCF files
        # for the patient.
        cell_vcf_paths = [(experimental_path / experimental_sample_id). \
                            with_suffix(".vcf") for experimental_sample_id 
                            in experimental_sample_ids]

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

        # TODO: Maybe remove this in Python 3.8.
        # This thread pool max-worker count is from the implementation in
        # Python 3.8. Assuming that Pathos adopts the same semantics,
        # this can be removed.
        with ThreadPool(min(32, os.cpu_count() + 4)) as pool:
            pool.map(process_cell, cell_vcf_paths)

    print("Running germline filter...")
    if processes > 1:
        with Pool(processes) as pool:
            list(tqdm(pool.imap(process_patient, all_germline_sample_ids),
                      total=len(all_germline_sample_ids), smoothing=0.01))
    else:
        list(map(process_patient, tqdm(all_germline_sample_ids, smoothing=0.1)))

    print("Done!")
