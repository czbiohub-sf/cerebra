import re
import numpy as np
from tqdm import tqdm

from ncls import NCLS

class GenomePosition():
	genome_pos_pattern = re.compile(r"(\d+):(\d+)-(\d+)")

	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.start = start
		self.end = end
	
	@classmethod
	def from_str(cls, pos_str):
		match = cls.genome_pos_pattern.match(pos_str)
		
		if not match:
			# FIXME: this is bad
			return cls("?", -1, -1)

		return cls(match[1], int(match[2]) - 1, int(match[3]))
	
	def __eq__(self, other):
		return self.chrom == other.chrom and self.start == other.start and self.end == other.end

	def __repr__(self):
		return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)

	def __str__(self):
		return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)
	
	def __len__(self):
		return self.end - self.start


class GenomeDataframeTree():
	def __init__(self, predicate, df):
		self.predicate = predicate
		self.df = df

		working_tree_map = {}

		# Iterating DataFrame rows :'(
		for idx, row in tqdm(df.iterrows(), total=len(df)):
			genome_pos = predicate(row)
			chrom = genome_pos.chrom

			if not chrom in working_tree_map:
				# (starts, ends, ids)
				working_tree_map[chrom] = ([], [], [])
			
			starts, ends, ids = working_tree_map[chrom]
			starts.append(genome_pos.start)
			ends.append(genome_pos.end)
			ids.append(idx)

		tree_map = {}

		for chrom, params in working_tree_map.items():
			starts, ends, ids = params
			tree_map[chrom] = NCLS(
				np.array(starts, dtype=np.long),
				np.array(ends, dtype=np.long),
				np.array(ids, dtype=np.long)
			)
		
		self.tree_map = tree_map
	
	def _compute_jaccard_index(self, pos_a, pos_b):
		range_a = range(pos_a.start - 1, pos_a.end)
		range_b = range(pos_b.start - 1, pos_b.end)

		if range_b.start > range_a.stop or range_a.start > range_b.stop:
			return 0
		
		intersection = range(max(range_a.start, range_b.start), min(range_a.stop, range_b.stop))

		# The following is equivalent to |A ∩ B| / |A ∪ B|, but avoids computing
		# a union.
		# |A ∩ B| / (|A| + |B| + |A ∩ B|)
		return len(intersection) / (len(range_a) + len(range_b) - len(intersection))
	
	def has_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return False
		
		# The overlap algorithm used by NCLS isn't inclusive of edges, so we pad
		# the edges of our query by 1.
		return tree.has_overlap(genome_pos.start - 1, genome_pos.end + 1)
	
	def get_first_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None
		
		starts = np.array([genome_pos.start], dtype=np.long)
		ends = np.array([genome_pos.end], dtype=np.long)
		ids = np.array([0], dtype=np.long)
		
		# The overlap algorithm used by NCLS isn't inclusive of edges, so we pad
		# the edges of our query by 1.
		_, row_ids = tree.first_overlap_both(starts - 1, ends + 1, ids)

		if len(row_ids) < 1:
			return None

		return self.df.iloc[row_ids[0]]

	def get_best_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None
		
		starts = np.array([genome_pos.start], dtype=np.long)
		ends = np.array([genome_pos.end], dtype=np.long)
		ids = np.array([0], dtype=np.long)
		
		# The overlap algorithm used by NCLS isn't inclusive of edges, so we pad
		# the edges of our query by 1.
		_, row_ids = tree.all_overlaps_both(starts - 1, ends + 1, ids)

		if len(row_ids) < 1:
			return None
		
		if len(row_ids) == 1:
			return self.df.iloc[row_ids[0]]
		
		rows = [self.df.iloc[row_id] for row_id in row_ids]
		scored_rows = [(row, self._compute_jaccard_index(genome_pos, self.predicate(row))) for row in rows]

		sorted_rows = sorted(scored_rows, key=lambda tup: tup[1], reverse=True)

		return sorted_rows[0][0]
	
	def get_first_containment(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None
		
		starts = np.array([genome_pos.start], dtype=np.long)
		ends = np.array([genome_pos.end], dtype=np.long)
		ids = np.array([0], dtype=np.long)
		
		_, row_ids = tree.all_containments_both(starts, ends, ids)

		if len(row_ids) < 1:
			return None

		return self.df.iloc[row_ids[0]]

def make_genome_pos_vcf(record):
	# Although including `chr` in the CHR column constitutes malformed VCF, it
	# may be present, so it should be removed.
	CHROM = record.CHROM.replace("chr", "")

	return GenomePosition(CHROM, record.affected_start, record.affected_end)

	# ref_len = len(record.REF)
	# alt_len = len(record.ALT)

	# if ref_len == 1 and alt_len == 1:
	# 	return GenomePosition(CHROM, POS, POS)
	# elif ref_len > 1 and alt_len == 1:
	# 	return GenomePosition(CHROM, POS, POS + ref_len)
	# elif alt_len > 1 and ref_len == 1:
	# 	return GenomePosition(CHROM, POS, POS + alt_len)
	# else: # multibase-for-multibase substitution
	# 	return GenomePosition(CHROM, POS, 1)

def make_genome_pos_gtf(record):
	return GenomePosition(record[0].replace("chr", ""), int(record[3]) - 1, int(record[4]))
