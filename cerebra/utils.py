import re

import numpy as np
import vcfpy
from hgvs import edit
from ncls import NCLS


class GenomePosition():
    genome_pos_pattern = re.compile(r"(.+):(\d+)-(\d+)")

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    @classmethod
    def from_str(cls, pos_str):
        match = cls.genome_pos_pattern.match(pos_str)

        if not match:
            return None

        return cls(match[1], int(match[2]) - 1, int(match[3]))

    @classmethod
    def from_vcf_record(cls, record):
        CHROM = record.CHROM.replace("chr", "")

        affected_ranges = [
            vcf_alt_affected_range(record.REF, alt) for alt in record.ALT
        ]
        start = record.begin + min(map(lambda r: r.start, affected_ranges),
                                   default=0)
        end = record.begin + max(map(lambda r: r.stop, affected_ranges),
                                 default=1)

        return cls(CHROM, start, end)

    @classmethod
    def from_vcf_record_pos(cls, record):
        CHROM = record.CHROM.replace("chr", "")
        return cls(CHROM, record.begin, record.begin + 1)

    @classmethod
    def from_gtf_record(cls, record):
        return cls(record[0].replace("chr", ""),
                   int(record[3]) - 1, int(record[4]))

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and \
        self.end == other.end

    def __repr__(self):
        return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)

    def __str__(self):
        return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)

    def __len__(self):
        return self.end - self.start

    def __contains__(self, other):
        same_chrom = other.chrom == self.chrom
        same_start = other.start >= self.start
        same_end = other.end <= self.end
        return same_chrom and same_start and same_end

    def __and__(self, other):
        if self.chrom != other.chrom:
            return None

        if other.start >= self.end or self.start >= other.end:
            return None

        return self.__class__(self.chrom, max(self.start, other.start),
                              min(self.end, other.end))

    # FIXME: This method may be too overloaded...
    def shifted_by(self, start, end=None):
        if isinstance(start, range):
            start, end = start.start, start.stop

        if end is None:
            end = start

        return self.__class__(self.chrom, self.start + start, self.end + end)

    def slice_within(self, other):
        if self.chrom != other.chrom:
            return None

        if self.start < other.start or self.end > other.end:
            return None

        return slice(self.start - other.start, self.end - other.start)


class GenomeIntervalTree():
    def __init__(self, predicate, records):
        self.predicate = predicate
        self.records = []

        working_tree_map = {}

        idx = 0
        for record in records:
            genome_pos = predicate(record)

            if genome_pos is None:
                continue

            chrom = genome_pos.chrom

            if chrom not in working_tree_map:
                # (starts, ends, ids)
                working_tree_map[chrom] = ([], [], [])

            starts, ends, ids = working_tree_map[chrom]
            starts.append(genome_pos.start)
            ends.append(genome_pos.end)
            ids.append(idx)

            self.records.append(record)
            idx += 1

        tree_map = {}

        for chrom, (starts, ends, ids) in working_tree_map.items():
            tree_map[chrom] = NCLS(np.array(starts, dtype=np.long),
                                   np.array(ends, dtype=np.long),
                                   np.array(ids, dtype=np.long))

        self.tree_map = tree_map

    def _intervals(self, chrom):
        return self.tree_map[chrom].intervals()

    def _make_query_params(self, genome_pos_list):
        starts = np.array([genome_pos.start for genome_pos in genome_pos_list])
        ends = np.array([genome_pos.end for genome_pos in genome_pos_list])
        ids = np.array(list(range(len(genome_pos_list))))

        return (starts, ends, ids)

    def _pick_best_record(self, from_ids=None, for_pos=None):
        if len(from_ids) < 1:
            return None

        if len(from_ids) == 1:
            return self.records[from_ids[0]]

        records = [self.records[record_id] for record_id in from_ids]

        scored_records = [(record,
                           self._compute_jaccard_index(for_pos,
                                                       self.predicate(record)))
                          for record in records]
        sorted_records = sorted(scored_records,
                                key=lambda tup: tup[1],
                                reverse=True)

        return sorted_records[0][0]

    def _compute_jaccard_index(self, pos_a, pos_b):
        intersection = pos_a & pos_b

        if not intersection:
            return 0

        # The following is equivalent to |A ∩ B| / |A ∪ B|, but avoids
        # computing a union.
        # |A ∩ B| / (|A| + |B| - |A ∩ B|)
        return len(intersection) / (len(pos_a) + len(pos_b) -
                                    len(intersection))

    def has_overlap(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return False

        return tree.has_overlap(genome_pos.start, genome_pos.end)

    def get_first_overlap(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return None

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.first_overlap_both(*qparams)

        if len(record_ids) < 1:
            return None

        return self.records[record_ids[0]]

    def get_best_overlap(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return None

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.all_overlaps_both(*qparams)

        return self._pick_best_record(from_ids=record_ids, for_pos=genome_pos)

    def get_all_overlaps(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return []

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.all_overlaps_both(*qparams)

        if any(map(lambda r: r >= len(self.records), record_ids)):
            print("uh oh i dont like this")
            print(len(self.records))
            print(record_ids)

        return [self.records[record_id] for record_id in record_ids]

    def get_first_containment(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return None

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.all_containments_both(*qparams)

        if len(record_ids) < 1:
            return None

        return self.records[record_ids[0]]

    def get_best_containment(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return None

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.all_containments_both(*qparams)

        return self._pick_best_record(from_ids=record_ids, for_pos=genome_pos)

    def get_all_containments(self, genome_pos):
        tree = self.tree_map.get(genome_pos.chrom)

        if not tree:
            return []

        qparams = self._make_query_params([genome_pos])
        _, record_ids = tree.all_containments_both(*qparams)

        return [self.records[record_id] for record_id in record_ids]


class GFFFeature():
    @classmethod
    def parse_gff_attributes(cls, attr_str):
        attr_dict = {}
        for key, value in (kv_str.split(' ')
                           for kv_str in re.split('; ?', attr_str) if kv_str):
            if '"' in value:
                value = value[1:-1]
            else:
                value = int(value)

            attr_dict[key] = value

        return attr_dict

    @property
    def is_forward_stranded(self):
        return self.strand == '+'

    @property
    def is_reverse_stranded(self):
        return self.strand == '-'

    def __init__(self, record):
        self.pos = GenomePosition.from_gtf_record(record)
        self.source = record[1]
        self.type = record[2]
        self.score = None if record[5] == '.' else float(record[5])
        self.strand = record[6]
        self.phase = None if record[7] == '.' else int(record[7])
        self.attributes = self.parse_gff_attributes(record[8])


def vcf_alt_affected_range(ref, alt):
    # TODO: This method currently only deals with simple substitutions.

    if alt.type in [vcfpy.SNV, vcfpy.MNV]:
        return range(len(ref))
    elif alt.type == vcfpy.INS:
        return range(2)
    elif alt.type == vcfpy.DEL:
        return range(1, len(ref))
    elif alt.type == vcfpy.INDEL:
        return range(len(ref))

    raise NotImplementedError()


def _seqs_are_equal(seq_a, seq_b, wildcard=None):
    if not len(seq_a) == len(seq_b):
        return False

    for a, b in zip(seq_a, seq_b):
        if a == wildcard or b == wildcard:
            continue

        if not a == b:
            return False

    return True


# This could be extended for other types of `SequenceVariant`s in the future if
# needed.
def sequence_variants_are_equivalent(seqvar_a,
                                     seqvar_b,
                                     strict_uncertain=False,
                                     strict_unknown=True,
                                     strict_silent=False):
    """Check if `seqvar_a` and `seqvar_b` are equivalent.
    Currently only works correctly for protein-level variants.

    Parameters
    ---------
    strict_uncertain : bool
        True if variant (position/edit) uncertainty is factored into
        this equivalency check. (default False)
    strict_unknown : bool
        True if unknown sequence units (e.g. 'X' for amino acids) should
        not match known sequence units. (default True)
    strict_silent : bool
        True if synonymous variants (e.g. 'Arg17=') should not match
        otherwise equivalent variants. (default False)
    """

    if not seqvar_a.ac == seqvar_b.ac:
        return False

    if not seqvar_a.type == seqvar_b.type:
        return False

    sv_type = seqvar_a.type

    if sv_type not in ["p"]:
        raise NotImplementedError()

    posedit_a, posedit_b = seqvar_a.posedit, seqvar_b.posedit

    if (posedit_a is None) or (posedit_b is None):
        return posedit_a is None and posedit_b is None

    if strict_uncertain and not posedit_a.uncertain == posedit_b.uncertain:
        return False

    pos_a, pos_b = posedit_a.pos, posedit_b.pos

    # TODO: Handle positional uncertainty

    if not pos_a == pos_b:
        return False

    edit_a, edit_b = posedit_a.edit, posedit_b.edit

    if not type(edit_a) is type(edit_b):
        print(type(edit_a), type(edit_b))
        return False

    _seqs_cmp = lambda a, b: _seqs_are_equal(
        a, b, wildcard=(None if strict_unknown else 'X'))

    if isinstance(edit_a, (edit.AARefAlt, edit.AAFs, edit.AAExt)):
        if (edit_a is None) or (edit_b is None):
            return edit_a is None and edit_b is None

        if not _seqs_cmp(edit_a.ref, edit_b.ref):
            return False

        if not _seqs_cmp(edit_a.alt, edit_b.alt):
            return False

        if strict_silent and (not edit_a.ref) and (not edit_a.alt):
            return False
    else:
        raise NotImplementedError()

    if isinstance(edit_a, (edit.AAFs, edit.AAExt)):
        if not edit_a.length == edit_b.length:
            return False

    if isinstance(edit_b, (edit.AAExt)):
        if not _seqs_cmp(edit_a.aaterm, edit_b.aaterm):
            return False

    return True
