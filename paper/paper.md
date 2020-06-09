---
title: 'cerebra: A tool for fast and accurate summarizing of variant calling format (.vcf) files'
tags:
    - python
    - genomics
    - variant calling
    - VCF
    - single-cell
    - cancer
authors:
 - name: Lincoln Harris
   orcid: 0000-0003-0872-7098
   affiliation: 1
 - name: Rohan Vanheusden
   orcid: 0000-0003-0872-7098
   affiliation: 1 
 - name: Olga Botvinnik
   orcid: 0000-0003-0872-7098
   affiliation: 1 
 - name: Spyros Darmanis
   orcid: 0000-0003-0872-7098
   affiliation: 1 
affiliations: 
- name: Chan Zuckerberg Biohub, San Francisco, CA
  index: 1
date: 20, August, 2020 
bibliography: paper.bib
---

## Motivation

A single "typo" in the genome can have massive consequences on an organism's biology.
Identifying the functional consequences of genomic typos (_i.e._ variants) is a fundamental challenge in bioinformatics. 
There exist tools for identifying variants and predicting their functional consequences, however, wrangling variant calls and functional predictions across thousands of samples represents an unsolved problem. 
`cerebra` addresses this need by offering a fast and accurate framework for summarizing variant calls and functional predictions across many samples. 

To find variants in the genome, researchers often begin with a [DNA-sequencing](https://en.wikipedia.org/wiki/DNA_sequencing) (DNA-seq) or [RNA-sequencing](https://en.wikipedia.org/wiki/RNA-Seq) (RNA-seq) experiment on their samples of interest.
After sequencing, the next step is alignment to the reference genome with tools like [STAR](https://github.com/alexdobin/STAR) or [BWA](http://bio-bwa.sourceforge.net/), followed by variant calling with tools like [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) 
or [freebayes](https://github.com/ekg/freebayes) [@star, @bwa, @haplocaller, @freebayes]. 
Variant callers produce tab delimited text files in the ([variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf), VCF)
for each processed sample, which encode the _genomic position_, _reference_ vs. _observed DNA sequence_, and _quality_
associated with each observed variant. 
An example VCF file:

```
##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM  POS ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
chr1	631391	.	C	T	72.28	.	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;QD=25.36;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:84,6,0
chr1	631862	.	G	A	286	.	AC=2;AF=1.00;AN=2;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=3.00;QD=28.73;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,8:8:24:300,24,0
chr1	1014274	rs8997	A	G	72.28	.	AC=2;AF=1.00;AN=2;DB;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;QD=30.97;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:84,6,0
chr1	1309460	.	G	A	245.98	.	AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=NaN;QD=27.24;SOR=4.174	GT:AD:DP:GQ:PL	1/1:0,7:7:21:260,21,0
```

Current methods for variant calling are incredibly powerful and robust, however, a single sequencing run can generate on the order of 10^8 unique VCF records, only a small portion of which are relevant to the researcher.

In addition, variant callers report only the genomic location and not the _functional_ consequences of the variant, _i.e._ the effect the variant has on the translated protein sequence.
We refer to these functional variants as "peptide-level variants." 
We introduce `cerebra`, a python package that provides fast and accurate peptide-level summarizing of VCF files.

## Functionality

`cerebra` comprises three modules: (1) `germline-filter` removes variants that are common between germline samples 
and samples of interest, (2) `count-mutations` reports total number of variants in each sample, and (3) `find-aa-mutations` reports likely peptide-level variants in each sample. 
Here we use _variant_ to refer to single nucleotide polymorphisms (SNPs) and short insertions and deletions. 
`cerebra` is not capable of reporting larger structural variants such as copy number variations and chromosomal rearrangements.

A data structure crucial to `cerebra` is the *genome interval tree*, which matches RNA transcripts
and peptides to each feature in the genome \autoref{fig1}. 
[Interval trees](https://en.wikipedia.org/wiki/Interval_tree) are self-balancing binary search trees that store numeric intervals and can quickly find every such interval that overlaps a given query interval (_[also see](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_). 
Given _n_ nodes, interval trees have theoretical average-case O(log*n*) and worst-case O(*n*) time complexity for search operations, making them tractable for genome-scale operations [@Cormen:2009, _[also see](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_].
Tree construction proceeds at O(*n*log*n*) time complexity, making construction rather than search the bottleneck for most VCF sets [@Alekseyenko:2007]. 
The _genome interval tree_ is constructed with a reference genome sequence ([FASTA format](https://en.wikipedia.org/wiki/FASTA_format), often with a `.fa` extension), and a genome annotation 
([gene transfer format, GTF](https://www.gencodegenes.org/pages/data_format.html) `.gtf` extension).
We rely on the [ncls](https://github.com/biocore-ntnu/ncls) library for fast interval tree construction and lookup operations.

We use [parallel processing](https://en.wikipedia.org/wiki/Multiprocessing) to stream in multiple VCF files at once. We extract relevant information -- including genomic interval, observed base, and read coverage -- from each variant record. In the `germline-filter` module variants are compared to one another and filtered out if found to be identical. In `count-mutations` variants are simply matched to whichever gene they came from. In `find-aa-mutations` variants are queried against our _genome interval tree_ -- if a matching interval is found we convert the DNA-level variant to a peptide-level variant. Eventually peptide-level variants from across all VCF are reported in tabular format. 

[todo: add in figure here]
