from collections import defaultdict
import re

import click
import numpy as np
import pandas as pd
import vcfpy
import hgvs.parser
from hgvs.sequencevariant import SequenceVariant
from tqdm import tqdm
from pathlib import Path
from pyfaidx import Fasta

from pathos.pools import _ProcessPool as Pool

from .protein_variant_predictor import ProteinVariantPredictor
from .utils import GenomePosition, GenomeIntervalTree, GFFFeature, vcf_alt_affected_range, sequence_variants_are_equal


class AminoAcidMutationFinder():
    def __init__(self, cosmic_df, annotation_df, genome_faidx):
        filtered_cosmic_df = self._make_filtered_cosmic_df(cosmic_df)

        print("Building COSMIC genome tree...")

        self._cosmic_genome_tree = GenomeIntervalTree(
            lambda row: GenomePosition.from_str(
                str(row["Mutation genome position"])),
            tqdm((record for idx, record in filtered_cosmic_df.iterrows()),
                 total=len(filtered_cosmic_df)))

        print("Building genome feature tree...")

        self._annotation_genome_tree = GenomeIntervalTree(
            lambda feat: feat.pos,
            tqdm((GFFFeature(row) for _, row in annotation_df.iterrows()),
                 total=len(annotation_df)))

        print("Building protein variant predictor...")

        self._protein_variant_predictor = ProteinVariantPredictor(
            self._annotation_genome_tree, genome_faidx)

    @classmethod
    def _make_filtered_cosmic_df(cls, cosmic_df):
        """Return a view of the input `cosmic_df` filtered for
        `Primary site` == `lung`."""

        filtered_db = cosmic_df.loc[cosmic_df["Primary site"] == "lung"]
        return filtered_db

    hgvs_parser = hgvs.parser.Parser(expose_all_rules=True)
    mutation_aa_delins_fix_pattern = re.compile(r"(p\.\d+_\d+)(\w+)>(\w+)")
    mutation_aa_uncertain_ins_del_fix_pattern = re.compile(
        r"(p\.)(\d+_\d+)(ins|del)(\d+)")

    @classmethod
    def _get_cosmic_record_protein_variant(cls, cosmic_record):
        mutation_aa = cosmic_record["Mutation AA"]
        accession = cosmic_record["Accession Number"]

        # Unfortunately, it seems that (some?) COSMIC mutation CDS strings are not quite HGVS-compliant; they use `>` rather than `del` + `ins`.
        mutation_aa = cls.mutation_aa_delins_fix_pattern.sub(
            r"\1del\2ins\3", mutation_aa)

        # COSMIC mutation CDS strings have uncertainty variants which are not HGVS-compliant.
        # mutation_cds = self.mutation_cds_uncertain_ins_del_fix_pattern.sub(r"\1(\2)\3(\4)", mutation_cds)

        cosmic_protein_posedit = (
            cls.hgvs_parser  # pylint: disable=no-member
            .parse_p_typed_posedit(mutation_aa)).posedit

        cosmic_protein_variant = SequenceVariant(
            ac=accession, type='p', posedit=cosmic_protein_posedit)

        return cosmic_protein_variant

    def _cosmic_subset_contains_genome_pos(self, genome_pos):
        return self._cosmic_genome_tree.has_overlap(genome_pos)

    def _make_mutation_counts_df(self, cell_genemuts_pairs):
        """Transform a list of tuples in the form
        [(cell_name, {gene_name: {aa_mutation}})] into a `DataFrame`."""

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
            aa_mutations = [list() for _ in range(len(cell_genemuts_pairs))]

            for idx, gene_muts in enumerate(cell_gene_muts):
                aa_mutations[idx] = list(gene_muts[gene_name])

            mutation_data[gene_name] = aa_mutations

        # TODO: cell_names index?
        return pd.DataFrame(index=cell_names, data=mutation_data)

    def find_cell_gene_aa_mutations(self, stream=None, path=None):
        """Create a `dict` mapping gene names to amino acid-level mutations
        found in that gene (filtered by COSMIC).
        Accepts a `stream` or a `path` as input."""

        vcf_reader = vcfpy.Reader.from_stream(
            stream) if stream is not None else vcfpy.Reader.from_path(path)

        gene_aa_mutations = defaultdict(set)

        for record in vcf_reader:
            record_pos = GenomePosition.from_vcf_record(record)

            # TODO: maybe clean this up into a method

            protein_variant_results = self._protein_variant_predictor \
                .predict_for_vcf_record(record)

            if not protein_variant_results:
                continue

            target_variants = []
            if self._cosmic_genome_tree:
                overlaps = self._cosmic_genome_tree.get_all_overlaps(
                    record_pos)

                target_variants = (
                    self._get_cosmic_record_protein_variant(overlap)
                    for overlap in overlaps)

                if not target_variants:
                    continue

            for result in protein_variant_results:
                predicted_variant = result.predicted_variant

                if target_variants:
                    for target_variant in target_variants:
                        if sequence_variants_are_equal(target_variant,
                                                       predicted_variant):
                            break
                    else:
                        # This protein variant didn't match any of the target
                        # variants; don't proceed.
                        continue

                # The predicted protein variant matches the protein variant
                # found in COSMIC.
                gene_name = result.transcript_feat.attributes["gene_name"]
                gene_aa_mutations[gene_name].add(str(predicted_variant))

                # The REFs are not the same; thus none of the ALTs can match.
                # Exit early.
                # if posedit.edit.ref != record.REF:
                #     continue

                # for alt in record.ALT:
                #     # TODO: Consider subset mutations in the future?

                #     affected_range = vcf_alt_affected_range(record.REF, alt)
                #     affected_pos = GenomePosition(
                #         record_pos.chrom,
                #         record_pos.start + affected_range.start,
                #         record_pos.end + affected_range.stop
                #     )

                #     # The mutations' positions are not the same.
                #     if overlap_pos != affected_pos:
                #         continue

                #     # The mutations' ALTs are not the same.
                #     if posedit.edit.alt != alt:
                #         continue

                #     # These mutations are the same!
                #     gene_aa_mutations.get(gene_name, d=[]).append(mutation_aa)

        return gene_aa_mutations

    def find_aa_mutations(self, paths, processes=4):
        """Create a `DataFrame` of mutation counts, where the row indices are
        cell names and the column indices are gene names."""
        def init_process(aa_mutation_finder):
            global current_process_aa_mutation_finder
            current_process_aa_mutation_finder = aa_mutation_finder

        def process_cell(path):
            return (
                Path(path).stem,
                current_process_aa_mutation_finder \
                    .find_cell_gene_aa_mutations(path=path)
            )

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
              default=1,
              prompt="number of processes to use for computation",
              type=int)
@click.option("--cosmicdb",
              prompt="path to cosmic db file (.tsv)",
              required=True)
@click.option("--annotation",
              prompt="path to genome annotation (.gtf)",
              required=True)
@click.option("--genomefa",
              prompt="path to full genome sequences (.fasta)",
              required=True)
@click.option("--outfile", prompt="path to output file (.csv)", required=True)
@click.argument("files", required=True, nargs=-1)
def find_aa_mutations(processes, cosmicdb, annotation, genomefa, outfile,
                      files):
    print("Beginning setup (this may take several minutes!)")

    print("Loading COSMIC database...")
    cosmic_df = pd.read_csv(cosmicdb, delimiter='\t')

    print("Loading genome annotation...")
    annotation_df = pd.read_csv(annotation, sep='\t', skiprows=5)

    print("Loading genome sequences...")
    genome_faidx = Fasta(genomefa)

    print("Building genome trees...")
    aa_mutation_finder = AminoAcidMutationFinder(cosmic_df, annotation_df,
                                                 genome_faidx)
    print("Setup complete.")

    print("Finding mutations...")
    result_df = aa_mutation_finder.find_aa_mutations(files,
                                                     processes=processes)

    print("Writing file...")
    outfile = Path(outfile)

    result_df.to_csv(outfile)

    print("Done!")
