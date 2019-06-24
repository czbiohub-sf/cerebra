import re
import numpy as np

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
			return None

		return cls(match[1], int(match[2]) - 1, int(match[3]))
	
	@classmethod
	def from_vcf_record(cls, record):
		# Although including `chr` in the CHR column constitutes malformed VCF, it
		# may be present, so it should be removed.
		CHROM = record.CHROM.replace("chr", "")

		return cls(CHROM, record.affected_start, record.affected_end)

	@classmethod
	def from_gtf_record(cls, record):
		return cls(record[0].replace("chr", ""), int(record[3]) - 1, int(record[4]))
	
	def contains(self, other):
		return other.chrom == self.chrom and other.start >= self.start and other.end <= self.end
	
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
		for idx, (_, row) in enumerate(df.iterrows()):
			genome_pos = predicate(row)

			if genome_pos is None:
				continue
			
			chrom = genome_pos.chrom

			if not chrom in working_tree_map:
				# (starts, ends, ids)
				working_tree_map[chrom] = ([], [], [])
			
			starts, ends, ids = working_tree_map[chrom]
			starts.append(genome_pos.start)
			ends.append(genome_pos.end)
			ids.append(idx)

		tree_map = {}

		for chrom, (starts, ends, ids) in working_tree_map.items():
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
		
		return tree.has_overlap(genome_pos.start, genome_pos.end)
	
	def get_first_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None
		
		starts = np.array([genome_pos.start], dtype=np.long)
		ends = np.array([genome_pos.end], dtype=np.long)
		ids = np.array([0], dtype=np.long)

		_, row_ids = tree.first_overlap_both(starts, ends, ids)

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
		
		_, row_ids = tree.all_overlaps_both(starts, ends, ids)

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
	
	def get_all_overlaps(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None
		
		starts = np.array([genome_pos.start], dtype=np.long)
		ends = np.array([genome_pos.end], dtype=np.long)
		ids = np.array([0], dtype=np.long)

		_, row_ids = tree.all_overlaps_both(starts, ends, ids)

		return [self.df.iloc[row_id] for row_id in row_ids]
