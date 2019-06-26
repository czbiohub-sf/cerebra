import click
import numpy as np
import pandas as pd
import vcfpy
import re
from tqdm import tqdm
from pathlib import Path

from pathos.pools import ProcessPool
from multiprocessing import Pool

from .utils import GenomePosition, GenomeIntervalTree

class MutationCounter():
    def __init__(self, cosmic_df, hg38_df):
        filtered_cosmic_df = self._make_filtered_cosmic_df(cosmic_df)
        
        self._cosmic_genome_tree = GenomeIntervalTree(
            lambda row: GenomePosition.from_str(str(row["Mutation genome position"])),
            (record for idx, record in filtered_cosmic_df.iterrows()))
        self._hg38_genome_tree = GenomeIntervalTree(
            GenomePosition.from_gtf_record,
            (record for idx, record in hg38_df.iterrows()))
    
    def _cosmic_subset_contains_genome_pos(self, genome_pos):
        return self._cosmic_genome_tree.has_overlap(genome_pos)

    def _find_containing_gene_record(self, genome_pos):
        return self._hg38_genome_tree.get_best_overlap(genome_pos)
    
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
        cell_gene_names = [list(gene_muts.keys()) for gene_muts in cell_gene_muts]
        # Flatten them into a set
        all_gene_names = set([gene_name for sublist in cell_gene_names for gene_name in sublist])
        
        mutation_data = {}
        for gene_name in all_gene_names:
            mutation_counts = np.zeros(len(cell_genemuts_pairs), dtype=int)

            for idx, gene_muts in enumerate(cell_gene_muts):
                mutation_counts[idx] = gene_muts.get(gene_name, 0)
            
            mutation_data[gene_name] = mutation_counts
        
        # TODO: cell_names index?
        return pd.DataFrame(index=cell_names, data=mutation_data)

    def find_cell_gene_mut_counts(self, stream=None, path=None):
        """Create a `dict` mapping gene names to the number of times a mutation
        was found in that gene. Accepts a `stream` or a `path` as input."""

        vcf_reader = vcfpy.Reader.from_stream(stream) if stream is not None else vcfpy.Reader.from_path(path)

        gene_mutation_counts = {}

        for record in vcf_reader:
            genome_pos = GenomePosition.from_vcf_record(record)
            # TODO: Filter out duplicates?
            # And is it better to filter out positional duplicates or gene name
            # duplicates?

            if not self._cosmic_subset_contains_genome_pos(genome_pos):
                continue
                        
            gene_row = self._find_containing_gene_record(genome_pos)

            if gene_row is None:
                continue

            gene_name = self._parse_gene_name(gene_row[8])

            gene_mutation_counts[gene_name] = gene_mutation_counts.get(gene_name, 0) + 1

        return gene_mutation_counts

    def find_mutation_counts(self, paths, processes=4):
        """Create a `DataFrame` of mutation counts, where the row indices are
        cell names and the column indices are gene names."""

        def process_cell(path):
            cell_name = Path(path).stem
            return (cell_name, self.find_cell_gene_mut_counts(path=path))
        
        cell_genemuts_pairs = None

        if processes > 1:
            with ProcessPool(processes) as pool:
                cell_genemuts_pairs = tqdm(pool.imap(process_cell, paths), total=len(paths), smoothing=0.1)
        else:
            cell_genemuts_pairs = map(process_cell, tqdm(paths))

        cell_genemuts_pairs = list(cell_genemuts_pairs)
        return self._make_mutation_counts_df(cell_genemuts_pairs)

@click.command()
@click.option("--processes", default=1, prompt="number of processes to use for computation", required=True, type=int)
@click.option("--cosmicdb", prompt="path to cosmic db file (.tsv)", required=True)
@click.option("--refgenome", prompt="path to reference genome (.gtf)", required=True)
@click.option("--outfile", prompt="path to output file (.csv)", required=True)
@click.argument("files", required=True, nargs=-1)
def count_mutations(processes, cosmicdb, refgenome, outfile, files):
    print("Beginning setup (this may take several minutes!)")

    print("Loading COSMIC database...")
    cosmic_df = pd.read_csv(cosmicdb, delimiter='\t')
    print("Loading reference genome...")
    refgenome_df = pd.read_csv(refgenome, delimiter='\t', header=None)
    
    print("Building genome trees...")
    mutation_counter = MutationCounter(cosmic_df, refgenome_df)
    print("Setup complete.")

    print("Finding mutations...")
    result_df = mutation_counter.find_mutation_counts(files, processes=processes)

    print("Writing file...")
    result_df.to_csv(outfile)
    
    print("Done!")
