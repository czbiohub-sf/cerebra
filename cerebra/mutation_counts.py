import numpy as np
import vcf
import os
import csv
import pandas as pd
import sys
import warnings
import click
import re
from tqdm import tqdm
from pathlib import Path

from pathos.pools import ProcessPool
from multiprocessing import Pool

from .utils import GenomePosition, GenomeDataframeTree

class MutationCounter():
    def __init__(self, cosmic_df, hg38_df):
        filtered_cosmic_df = self._make_filtered_cosmic_df(cosmic_df)
        
        print("Building COSMIC genome tree...")
        self._cosmic_genome_tree = GenomeDataframeTree(lambda row: GenomePosition.from_str(str(row["Mutation genome position"])), filtered_cosmic_df)
        print("Building reference genome tree...")
        self._hg38_genome_tree = GenomeDataframeTree(GenomePosition.from_gtf_record, hg38_df)
        print("Setup complete.")
    
    def _cosmic_subset_contains_genome_pos(self, genome_pos):
        return self._cosmic_genome_tree.has_overlap(genome_pos)

    def _find_containing_gene_record(self, genome_pos):
        return self._hg38_genome_tree.get_best_overlap(genome_pos)
    
    gene_name_pattern = re.compile(r"gene_name \"(.+?)\"")
    def _parse_gene_name(self, record):
        gene_name_match = self.gene_name_pattern.search(record[8])

        if not gene_name_match:
            return None

        return gene_name_match[1]

    def _make_filtered_cosmic_df(self, cosmic_df):
        """ returns the COSMIC database after lung adeno filter """
        print("Filtering COSMIC database...")
        filtered_by_primary_site = cosmic_df.index[cosmic_df["Primary site"] == "lung"].tolist()
        database_filter = cosmic_df.iloc[filtered_by_primary_site]

        return database_filter

    def _make_mutation_counts_df(self, cell_genemuts_pairs):
        """ creates the cell/mutation counts table from the raw output that 
            get_gene_cell_muts_counts provides """

        cell_names, cell_gene_muts = list(zip(*cell_genemuts_pairs))
        cell_names, cell_gene_muts = list(cell_names), list(cell_gene_muts)

        # Aggregate all of the gene names into a 2D array
        cell_gene_names = [list(gene_muts.keys()) for gene_muts in cell_gene_muts]
        # Flatten them into a set
        all_gene_names = set([gene_name for sublist in cell_gene_names for gene_name in sublist])
        
        mutation_data = {}
        for gene_name in all_gene_names:
            mutation_counts = np.zeros(len(cell_genemuts_pairs), dtype=np.int32)

            for idx, gene_muts in enumerate(cell_gene_muts):
                mutation_counts[idx] = gene_muts.get(gene_name, 0)
            
            mutation_data[gene_name] = mutation_counts
        
        # TODO: cell_names index?
        return pd.DataFrame(index=cell_names, data=mutation_data)

    def find_cell_gene_mut_counts(self, filename):
        """ returns a dictionary mapping a gene name to the number of observed
            mutations """

        # PyVCF documentation claims that it automatically infers compression type
        # from the file extension.
        vcf_reader = vcf.Reader(filename=filename)

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
                #print("gene_name miss")
                continue

            gene_name = self._parse_gene_name(gene_row)

            # gene_name = find_gene_name(genome_pos, hg38_genome_tree)

            gene_mutation_counts[gene_name] = gene_mutation_counts.get(gene_name, 0) + 1

        return gene_mutation_counts

    def find_mutation_counts(self, filenames, threads=4):
        """ generate a cell x gene mutation counts table from a set of germline filtered vcfs """

        def process_cell(filename):
            cell_name = Path(filename).stem
            return (cell_name, self.find_cell_gene_mut_counts(filename))
        
        cell_genemuts_pairs = None

        if threads > 1:
            with ProcessPool(threads) as pool:
                cell_genemuts_pairs = tqdm(pool.imap(process_cell, filenames), total=len(filenames), smoothing=0.1)
        else:
            cell_genemuts_pairs = map(process_cell, tqdm(filenames))

        cell_genemuts_pairs = list(cell_genemuts_pairs)
        return self._make_mutation_counts_df(cell_genemuts_pairs)

@click.command()
@click.option("--threads", default=4, prompt="number of threads", required=True, type=int)
@click.option("--cosmicdb", prompt="path to cosmic db file (.tsv)", required=True)
@click.option("--refgenome", prompt="path to reference genome (.gtf)", required=True)
@click.option("--outfile", prompt="path to output file (.csv)", required=True)
@click.argument("files", required=True, nargs=-1)
def count_mutations(threads, cosmicdb, refgenome, outfile, files):
    print("Beginning setup (this may take several minutes!)")

    print("Loading COSMIC database...")
    cosmic_df = pd.read_csv(cosmicdb, delimiter='\t')
    print("Loading reference genome...")
    refgenome_df = pd.read_csv(refgenome, delimiter='\t', header=None)
    
    mutation_counter = MutationCounter(cosmic_df, refgenome_df)

    print("Finding mutations...")
    result_df = mutation_counter.find_mutation_counts(files, threads=threads)

    print("Writing file...")
    result_df.to_csv(outfile)
    
    print("Done!")
