from collections import defaultdict
from pathlib import Path

import click
import numpy as np
import pandas as pd
import vcfpy
import re
from pathos.pools import _ProcessPool as Pool
from tqdm import tqdm

from .utils import GenomePosition, GenomeIntervalTree


class MutationCounter():
    def __init__(self, cosmic_df, annotation_df):
        if cosmic_df is not None:
            filtered_cosmic_df = self._make_filtered_cosmic_df(cosmic_df)

            self._cosmic_genome_tree = GenomeIntervalTree(
                lambda row: GenomePosition.from_str(
                    str(row["Mutation genome position"])),
                (row for idx, row in filtered_cosmic_df.iterrows()))
        else:
            self._cosmic_genome_tree = None

        self._annotation_genome_tree = GenomeIntervalTree(
            lambda feat: GenomePosition.from_gtf_record(feat),
            (row for _, row in annotation_df.iterrows()))

    def _filter_includes_genome_pos(self, genome_pos):
        if not self._cosmic_genome_tree:
            return True

        return self._cosmic_genome_tree.has_overlap(genome_pos)

    def _find_containing_gene_record(self, genome_pos):
        return self._annotation_genome_tree.get_best_overlap(genome_pos)

    gene_name_pattern = re.compile(r"gene_name \"(.+?)\"")

    def _parse_gene_name(self, metadata):
        """Parse out a gene name from a GTF metadata string."""

        gene_name_match = self.gene_name_pattern.search(metadata)

        if not gene_name_match:
            return None

        return gene_name_match[1]

    def _make_filtered_cosmic_df(self, cosmic_df):
        """Return a view of the input `cosmic_df` filtered for
        `Primary site` == `lung`."""

        filtered_db = cosmic_df.loc[cosmic_df["Primary site"] == "lung"]
        return filtered_db

    def _make_mutation_counts_df(self, cell_genemuts_pairs):
        """Transform a list of tuples in the form
        (cell_name, {gene_name: mutation_count}) into a `DataFrame`."""

        cell_names, cell_gene_muts = list(zip(*cell_genemuts_pairs))
        cell_names, cell_gene_muts = list(cell_names), list(cell_gene_muts)

        # Aggregate all of the gene names into a 2D array
        cell_gene_names = [
            list(gene_muts.keys()) for gene_muts in cell_gene_muts
        ]

        # Flatten them into a set
        all_gene_names = set([
            gene_name for sublist in cell_gene_names for gene_name in sublist
        ])

        mutation_data = {}
        for gene_name in all_gene_names:
            mutation_counts = np.zeros(len(cell_genemuts_pairs), dtype=int)

            for idx, gene_muts in enumerate(cell_gene_muts):
                mutation_counts[idx] = gene_muts[gene_name]

            mutation_data[gene_name] = mutation_counts

        # TODO: cell_names index?
        return pd.DataFrame(index=cell_names, data=mutation_data)

    def find_cell_gene_mut_counts(self, stream=None, path=None):
        """Create a `dict` mapping gene names to the number of times a mutation
        was found in that gene. Accepts a `stream` or a `path` as input."""

        vcf_reader = vcfpy.Reader.from_stream(
            stream) if stream is not None else vcfpy.Reader.from_path(path)

        gene_mutation_counts = defaultdict(int)

        for record in vcf_reader:
            genome_pos = GenomePosition.from_vcf_record(record)
            # TODO: Filter out duplicates?
            # And is it better to filter out positional duplicates or gene name
            # duplicates?

            # TODO: Report all relevant genes, not just "best" one.
            gene_row = self._find_containing_gene_record(genome_pos)

            if gene_row is None:
                continue

            gene_name = self._parse_gene_name(gene_row[8])

            if gene_name is None:
                continue

            if self._filter_includes_genome_pos(genome_pos):
                continue

            gene_mutation_counts[
                gene_name] = gene_mutation_counts[gene_name] + 1

        return gene_mutation_counts

    def find_mutation_counts(self, paths, processes=1):
        """Create a `DataFrame` of mutation counts, where the row indices are
        cell names and the column indices are gene names."""
        def init_process(mutation_counter):
            global current_process_mutation_counter
            current_process_mutation_counter = mutation_counter

        def process_cell(path):
            cell_name = Path(path).stem
            return (cell_name,
                    current_process_mutation_counter.find_cell_gene_mut_counts(
                        path=path))

        if processes > 1:
            with Pool(processes, initializer=init_process,
                      initargs=(self, )) as pool:
                results = list(
                    tqdm(pool.imap(process_cell, paths),
                         total=len(paths),
                         smoothing=0.01))
        else:
            init_process(self)
            results = list(map(process_cell, tqdm(paths)))

        return self._make_mutation_counts_df(results)


@click.command()
@click.option("--processes",
              "num_processes",
              help="number of processes to use for computation",
              default=1)
@click.option("--cosmicdb",
              "cosmicdb_path",
              help="optional path to cosmic db file (.tsv)",
              default=None)
@click.option("--annotation",
              "annotation_path",
              help="path to reference genome (.gtf)",
              required=True)
@click.option("--output",
              "output_path",
              help="path to output file (.csv)",
              required=True)
@click.argument("input_files", required=True, nargs=-1)
def count_mutations(num_processes, cosmicdb_path, annotation_path, output_path,
                    input_files):
    print("Beginning setup (this may take several minutes!)")

    if cosmicdb_path:
        print("Loading COSMIC database...")
        cosmic_df = pd.read_csv(cosmicdb_path, sep='\t')
    else:
        cosmic_df = None

    print("Loading genome annotation...")
    annotation_df = pd.read_csv(annotation_path, sep='\t')

    print("Building genome trees...")
    mutation_counter = MutationCounter(cosmic_df, annotation_df)
    print("Setup complete.")

    print("Finding mutations...")
    result_df = mutation_counter.find_mutation_counts(input_files,
                                                      processes=num_processes)

    print("Writing file...")
    output_path = Path(output_path)

    result_df.to_csv(output_path)

    print("Done!")
