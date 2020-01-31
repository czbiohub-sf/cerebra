from collections import defaultdict, namedtuple
from itertools import tee
import re 

from Bio import Alphabet
from Bio.Seq import Seq  # need to clean this up
from hgvs.utils.altseq_to_hgvsp import AltSeqToHgvsp

from .utils import GenomePosition, GenomeIntervalTree, vcf_alt_affected_range

_VariantTreeTranscriptRecord = namedtuple("VariantTreeTranscriptRecord", [
    "feat", "five_prime_cds_feat", "three_prime_cds_feat",
    "five_prime_utr_feats", "three_prime_utr_feats", "translated_feats"
])

_VariantTreeCodingRecord = namedtuple("VariantTreeCodingRecord",
                                      ["transcript", "feat"])

_LocalizedSeq = namedtuple("LocalizedSeq", ["seq", "pos"])

# Create a mock `RefTranscriptData` type which behaves sufficiently similarly
# for the purposes of this module. Used to avoid the initialization step of the
# actual type.
_MockRefTranscriptData = namedtuple("RefTranscriptData", [
    "transcript_sequence", "aa_sequence", "cds_start", "cds_stop",
    "protein_accession"
])

# Create a mock `AltTranscriptData` type which behaves sufficiently similarly
# for the purposes of this module. Used to avoid the initialization step of the
# actual type.
_MockAltTranscriptData = namedtuple("AltTranscriptData", [
    "transcript_sequence", "aa_sequence", "cds_start", "cds_stop",
    "protein_accession", "is_frameshift", "variant_start_aa",
    "frameshift_start", "is_substitution", "is_ambiguous"
])

ProteinVariantResult = namedtuple(
    "ProteinVariantResult",
    ["query_variant", "predicted_variant", "transcript_feat"])


class ProteinVariantPredictor():
    # TODO: Now I think it would be nice if this was written to use GFF3
    # instead of GTF because of the hierarchic feature structure GFF3 provides,
    # which is is already being emulated (to an extent) below.
    def __init__(self, annotation_genome_tree, genome_faidx):
        self.genome_fasta = genome_faidx

        transcript_feats_dict = defaultdict(lambda: defaultdict(list))
        for feature in annotation_genome_tree.records:
            if feature.attributes.get("transcript_type") != "protein_coding":
                continue

            feature_type = feature.type
            if feature_type not in [
                "transcript", "CDS", "UTR", "start_codon", "stop_codon"
            ]:
                continue

            transcript_id = feature.attributes["transcript_id"]
            transcript_feats_dict[transcript_id][feature_type].append(feature)

        transcript_records = {}
        cds_records = []
        for transcript_id, features in transcript_feats_dict.items():
            # There should only be one transcript record...
            transcript_feat = features["transcript"][0]

            # Hopefully redundant check
            if not transcript_feat:
                print(f"Unexpected missing transcript feature for transcript"
                      " ID {transcript_id}")
                continue

            cds_feats = features["CDS"]

            # Sort the CDS features so that the CDS adjacent to the 5' UTR is
            # first and the CDS adjacent to the 3' UTR is last.
            is_forward_stranded = transcript_feat.is_forward_stranded
            cds_feats.sort(reverse=(not is_forward_stranded),
                           key=lambda feat: feat.pos.start
                           if is_forward_stranded else feat.pos.end)

            five_prime_cds_feat = cds_feats[0]
            three_prime_cds_feat = cds_feats[-1]

            utr_feats = features["UTR"]

            is_before = lambda rh: lambda lh: lh.pos.end <= rh.pos.start
            is_after = lambda rh: lambda lh: is_before(lh)(rh)

            utr_filter_5p = is_before(
                five_prime_cds_feat) if is_forward_stranded else is_after(
                five_prime_cds_feat)
            utr_filter_3p = is_after(
                three_prime_cds_feat) if is_forward_stranded else is_before(
                three_prime_cds_feat)

            five_prime_utr_feats = filter(utr_filter_5p, utr_feats)
            three_prime_utr_feats = filter(utr_filter_3p, utr_feats)

            transcript_record = _VariantTreeTranscriptRecord(
                feat=transcript_feat,
                five_prime_cds_feat=five_prime_cds_feat,
                three_prime_cds_feat=three_prime_cds_feat,
                five_prime_utr_feats=five_prime_utr_feats,
                three_prime_utr_feats=three_prime_utr_feats,
                translated_feats=features["start_codon"] + cds_feats +
                                 features["stop_codon"])

            transcript_records[transcript_id] = transcript_record

            for cds_feat in cds_feats:
                cds_records.append(
                    _VariantTreeCodingRecord(transcript=transcript_record,
                                             feat=cds_feat))

        self.transcript_records = transcript_records
        self.tree = GenomeIntervalTree(lambda record: record.feat.pos,
                                       cds_records)

        if 'chrM' in self.genome_fasta.keys():
            self.genome_fasta_mito = genome_faidx["chrM"][:]
        else:
            self.genome_fasta_mito = None

    @classmethod
    def _merged_and_sorted_intersecting_intervals(cls, intervals):
        current_interval = None

        for interval in sorted(intervals):
            if current_interval is None:
                current_interval = interval
                continue

            i0, i1, c0, c1 = *interval, *current_interval

            if i0 > c1:  # i1 < c0 shouldn't be possible due to sort
                # This new interval doesn't overlap with the current one; yield
                # the old interval and rotate.
                yield current_interval
                current_interval = interval
                continue

            # This interval intersects with the current one; merge them
            current_interval = (min(i0, c0), max(i1, c1))

        if current_interval is not None:
            yield current_interval

    @classmethod
    def _splice_seq(cls, seq, intervals):
        ordered_intervals = cls._merged_and_sorted_intersecting_intervals(
            intervals)

        # Create empty `Seq` with same alphabet
        spliced_seq = seq[0:0]

        for interval in ordered_intervals:
            spliced_seq += seq[slice(*interval)]

        return spliced_seq

    def predict_for_vcf_record(self, vcf_record):

        def check_str(s):
            ''' check a string to see if it has non-sequence characters '''
            match = re.match("^[AGCTN]*$", s)
            return match is not None


        def handle_mito(_tx_pos):
            ''' handles this chrM wierdness '''
            mito_seq = None
            if self.genome_fasta_mito:
                mito_seq = Seq(self.genome_fasta_mito.seq[_tx_pos.start:_tx_pos.end], 
                            alphabet=Alphabet.generic_dna)
            return mito_seq


        variant_results = []
        record_pos = GenomePosition.from_vcf_record_pos(vcf_record)

        ref = vcf_record.REF
        for alt in vcf_record.ALT:
            # Create a GenomePosition representing the range affected by the
            # ALT sequence.
            affected_pos = record_pos.shifted_by(
                vcf_alt_affected_range(ref, alt))

            coding_overlaps = self.tree.get_all_overlaps(affected_pos)
            transcript_ids = set(
                record.transcript.feat.attributes["transcript_id"]
                for record in coding_overlaps)

            for tx_id in transcript_ids:
                transcript = self.transcript_records[tx_id]
                tx_pos = transcript.feat.pos
                ref_pos = record_pos.shifted_by(0, len(ref) - 1)

                # We want to make sure that the transcript sequence we pull
                # includes the REF so that we can cleanly replace it with the
                # ALT later.
                tx_pos = tx_pos.shifted_by(
                    min(0, ref_pos.start - tx_pos.start),
                    max(0, ref_pos.end - tx_pos.end))

                ref_tx_seq = Seq(self.genome_fasta["chr" + tx_pos.chrom]
                                 [tx_pos.start:tx_pos.end].seq,
                                 alphabet=Alphabet.generic_dna)

                if not check_str(str(ref_tx_seq)): # non-DNA seq chars
                    if tx_pos.chrom == 'M':
                        ref_tx_seq = handle_mito(tx_pos)
                        if not ref_tx_seq:
                            continue
                    else:   # this shouldn't happen; skip iter
                        continue 

                ref_slice = ref_pos.slice_within(tx_pos)

                # TODO: The '*' ALT indicates that the allele was missing.
                # Right now, this is treated as a deletion of the entire REF
                # sequence, but perhaps this should be handled differently...
                alt_value = '' if alt.value == '*' else alt.value

                alt_tx_seq = ref_tx_seq[:ref_slice.start] \
                             + alt_value \
                             + ref_tx_seq[ref_slice.stop:]

                alt_tx_seq_len_delta = len(alt_tx_seq) - len(ref_tx_seq)

                ref_splice_intervals, alt_splice_intervals = tee(
                    (feat.pos.slice_within(tx_pos).indices(len(tx_pos))[:2]
                     for feat in transcript.translated_feats), 2)

                record_pos_tx_slice = record_pos.slice_within(tx_pos)

                # NOTE: Because of the way the splicing interval offsets are
                # calculated, insertions should not affect the length of the
                # ALT coding sequence; e.g. an insertion which starts before
                # the 5' CDS will just shift the splicing intervals over by its
                # length.
                alt_splice_intervals = (
                    # If the variant comes after a feature's position, the
                    # interval doesn't need to be changed; if it does,
                    # the length difference between the ALT and REF is added.
                    # This applies individually for both the start and end of
                    # each interval.
                    (ival[0] + (0 if record_pos_tx_slice.start >= ival[0] else
                                alt_tx_seq_len_delta),
                     ival[1] + (0 if record_pos_tx_slice.start >= ival[1] else
                                alt_tx_seq_len_delta))
                    for ival in alt_splice_intervals)

                ref_splice_intervals = list(ref_splice_intervals)
                alt_splice_intervals = list(alt_splice_intervals)

                ref_coding_seq = self._splice_seq(ref_tx_seq,
                                                  ref_splice_intervals)
                alt_coding_seq = self._splice_seq(alt_tx_seq,
                                                  alt_splice_intervals)

                if transcript.feat.is_reverse_stranded:
                    ref_coding_seq = ref_coding_seq.reverse_complement()
                    alt_coding_seq = alt_coding_seq.reverse_complement()

                # TODO: For instances where start/stop codons have been
                # altered, consider splicing in UTRs up until start/stop
                # codons.

                n_term_padding = -(transcript.five_prime_cds_feat.phase % -3)

                ref_c_term_trim = (n_term_padding + len(ref_coding_seq)) % 3
                alt_c_term_trim = (n_term_padding + len(alt_coding_seq)) % 3
                ref_trim_slice = slice(0,
                                       len(ref_coding_seq) - ref_c_term_trim)
                alt_trim_slice = slice(0,
                                       len(alt_coding_seq) - alt_c_term_trim)

                ref_coding_seq = (
                                     'N' * n_term_padding) + ref_coding_seq[
                                     ref_trim_slice]
                alt_coding_seq = (
                                     'N' * n_term_padding) + alt_coding_seq[
                                     alt_trim_slice]

                ref_aa_seq = ref_coding_seq.translate()
                alt_aa_seq = alt_coding_seq.translate()

                protein_id = transcript.feat.attributes["protein_id"]

                # Alt transcript data flags

                is_frameshift = (len(ref_coding_seq) -
                                 len(alt_coding_seq)) % 3 != 0
                variant_start_aa = None
                is_substitution = False
                # The implementation for `is_ambiguous` used by `hgvs`
                # considers "ambiguity" to mean "more than one stop codon."
                is_ambiguous = alt_aa_seq.count('*') > 1

                min_coding_seq_length = min(len(ref_coding_seq),
                                            len(alt_coding_seq))
                for pos in range(min_coding_seq_length):
                    ref_na, alt_na = ref_coding_seq[pos], alt_coding_seq[pos]

                    if ref_na != alt_na:
                        variant_start_aa = pos // 3
                        break
                else:
                    if len(ref_coding_seq) == len(alt_coding_seq):
                        # The coding sequence was unchanged by the variant.
                        continue

                    variant_start_aa = min_coding_seq_length // 3

                # NOTE: Not sure if `hgvs` expects `is_substitution` to be
                # evaluated on the amino acid-level or the nucleic acid-level.
                # Currently evaluated on the amino acid-level; therefore silent
                # substitutions aren't taken into account.
                min_aa_seq_length = min(len(ref_aa_seq), len(alt_aa_seq))
                if len(ref_aa_seq) == len(alt_aa_seq):
                    for pos in range(min_aa_seq_length):
                        ref_aa, alt_aa = ref_aa_seq[pos], alt_aa_seq[pos]

                        if ref_aa != alt_aa:
                            if is_substitution:
                                is_substitution = False
                                break

                            is_substitution = True

                if not is_substitution and (
                    (len(ref_aa_seq) > len(alt_aa_seq)
                     and ref_aa_seq.startswith(alt_aa_seq)) or
                    (len(alt_aa_seq) > len(ref_aa_seq)
                     and alt_aa_seq.startswith(ref_aa_seq))):
                    # FIXME (by filing a bug report?):
                    # `hgvs` has a weird bug where deletions/insertions that do
                    # not cause any observable variation at the protein level
                    # cause an `IndexError` while attempting to find the first
                    # varying amino acid. These must be skipped for now.
                    continue

                ref_tx_data = _MockRefTranscriptData(
                    transcript_sequence=str(ref_coding_seq),
                    aa_sequence=str(ref_aa_seq),
                    cds_start=None,  # Only used for initialization
                    cds_stop=None,  # Only used for initialization
                    protein_accession=protein_id)
                alt_tx_data = _MockAltTranscriptData(
                    transcript_sequence=str(alt_coding_seq),
                    aa_sequence=str(alt_aa_seq),
                    cds_start=None,  # Only used for initialization
                    cds_stop=None,  # Only used for initialization
                    protein_accession=protein_id,
                    is_frameshift=is_frameshift,
                    variant_start_aa=variant_start_aa + 1,  # 1-indexed
                    frameshift_start=None,  # Seems to be unused currently
                    is_substitution=is_substitution,
                    # The `AltTranscriptData` builder used by `hgvs` sets
                    # `is_ambiguous` to true if there are multiple stop codons
                    # in the AA sequence. This doesn't work in this code's
                    # context, because frameshifts typically create many new
                    # stop codons. It is uncertain why `hgvs` doesn't run into
                    # this issue; perhaps it's because it indiscriminately pads
                    # sequences that can't be cleanly translated.
                    is_ambiguous=(not is_frameshift) and is_ambiguous)

                protein_variant_builder = AltSeqToHgvsp(
                    ref_tx_data, alt_tx_data)

                protein_variant = protein_variant_builder.build_hgvsp()

                variant_results.append(
                    ProteinVariantResult(query_variant=alt,
                                         predicted_variant=protein_variant,
                                         transcript_feat=transcript.feat))

        return variant_results
