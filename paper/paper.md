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
or [freebayes](https://github.com/ekg/freebayes) [@star; @bwa; @haplocaller; @freebayes]. 
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
and peptides to each feature in the genome (\autoref{workflow}). 
[Interval trees](https://en.wikipedia.org/wiki/Interval_tree) are self-balancing binary search trees that store numeric intervals and can quickly find every such interval that overlaps a given query interval (_[also see](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_). 
Given _n_ nodes, interval trees have theoretical average-case O(log*n*) and worst-case O(*n*) time complexity for search operations, making them tractable for genome-scale operations [@Cormen:2009, _[also see](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_].
Tree construction proceeds at O(*n*log*n*) time complexity, making construction rather than search the bottleneck for most VCF sets [@Alekseyenko:2007]. 
The _genome interval tree_ is constructed with a reference genome sequence ([FASTA format](https://en.wikipedia.org/wiki/FASTA_format), often with a `.fa` extension), and a genome annotation 
([gene transfer format, GTF](https://www.gencodegenes.org/pages/data_format.html) `.gtf` extension).
We rely on the [ncls](https://github.com/biocore-ntnu/ncls) library for fast interval tree construction and lookup operations.

We use [parallel processing](https://en.wikipedia.org/wiki/Multiprocessing) to stream in multiple VCF files at once. We extract relevant information -- including genomic interval, observed base, and read coverage -- from each variant record. In the `germline-filter` module variants are compared to one another and filtered out if found to be identical. In `count-mutations` variants are simply matched to whichever gene they came from. In `find-aa-mutations` variants are queried against our _genome interval tree_ -- if a matching interval is found we convert the DNA-level variant to a peptide-level variant. Eventually peptide-level variants from across all VCF are reported in tabular format. 

![Workflow describing the `find-aa-mutations` module. We construct a genome interval tree from a genome annotation (.gtf) and a reference genome sequence (.fa), then processing VCF files in parallel to create a single tabular output file.\label{workflow}](fig1.jpg)

## `germline-filter`

If the research project is centered around a "tumor/pathogenic vs control" question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between the control and the experimental tissue so as to not bias the results by including non-pathogenic variants. 
The user provides a very simple metadata file that indicates which experimental samples correspond to which control samples.
Using the [vcfpy](https://pypi.org/project/vcfpy/) library we quickly identify shared variants across control/experimental matched VCF files, then write new VCFs that contain only the unique variants. 
These steps are performed by a [subprocess pool](https://pypi.org/project/pathos/) so that we can quickly process "chunks" of input in a [parallel](https://en.wikipedia.org/wiki/Multiprocessing) manner. 
There is also the option to limit the reported variants to those found in NCBI's [dbSNP](https://www.ncbi.nlm.nih.gov/books/NBK21088/) and the Wellcome Sanger Institute's [COSMIC](https://cancer.sanger.ac.uk/cosmic) databases. 
This option is designed to give the user a higher degree of confidence in the pathogenic nature of each variant -- if independent experiments have reported a given variant in human tissue, there is a higher likilihood that it is pathogenic. 
The output of `germline-filter` is a set of trimmed-down VCF files. 

If you have access to "control" tissue and your experimental question is concerned with differences between tumor/pathogenic tissue and control tissue, then `germline-filter` is the right place to start.
`germline-filter` will produce a new set of VCFs, which you'll use for the next two steps.
If you do not have access to "control" tissue, then proceed directly to `count-mutations` or `find-aa-mutations`.

## `count-mutations`
The `count-mutations` module reports the raw variant counts for every gene across every sample.
We first create a _genome interval tree_ from the reference GTF, then read in a VCF file and convert it to a [vcfpy](https://pypi.org/project/vcfpy/) object, then processes VCF records in [parallel](https://en.wikipedia.org/wiki/Multiprocessing). 
Each variant is matched to its corresponding gene, and gene-wise counts are stored in shared memory. 
We then report the raw number of variants found in each sample. 
The output is a CSV file that contains counts for each sample versus every gene in the genome. 

## `find-aa-mutations`
The `find-aa-mutations` module reports the peptide-level consequence of variants in the genome.
First we load the reference GTF, then construct an index (.fai) of the genome fasta file with [pyfaidx](https://pypi.org/project/pyfaidx/) to enable fast random memory access. 
We then create a _genome interval tree_ that will be used to quickly match genomic coordinates from VCF records to peptide-level variants. 
If working  with cancer samples, the user has the option to filter out all variants that are not found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database and are therefore unlikely to be pathogenic.

VCF files are read in simultaneously; individual records are converted to _GenomePosition_ objects to keep track of their genomic intervals and observed DNA bases.
_GenomePositions_ are then queried against the _genome interval tree_. 
If an overlapping interval is found we retrieve the peptide-level variant from this node of the _genome interval tree_. 
Peptide-level variants are converted to [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, 
in acordance to the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature. 
The output is a heirarchically ordered text file (CSV or JSON) that reports the the Ensemble protein ID and the gene associated with each variant, for each experimental sample. 

We should stress that `find-aa-mutations` does not *definitively* report peptide-level variants but rather the *likely*
set of peptide variants. 
Definitively reporting protein variants requires knowledge of alternate splicing -- this represents an open problem in scRNA-seq [@Huang:2017]. 
For example, if a read picks up a variant in exon 2 of geneA, we can report each of the potential spliceforms of geneA that contain exon 2, but we **cannot** infer which of those particular spliceforms are actually present in our sample. 
Thus we report all possible spliceforms; determining the spliceform landscape of an individual cell from scRNA-seq is outside the scope of this project. 

We tested `find-aa-mutations` on a set of high-quality reference-grade VCF files from the [Genome in a Bottle consortium](https://www.nist.gov/programs-projects/genome-bottle) [@GiaB_orig; @GiaB_adnl]. 
Each of the seven VCF files was quite large, (~2GB) and `cerebra` was run on standard hardware (MacBook Pro, 2.5GHz quad-core processor, 16 GB RAM). 
`cerebra` processed the seven files in 44 minutes, see *Figure 2*. 
The Genome in a Bottle VCFs are much larger than the VCFs generated by a typical sequencing experiment; we thought it prudent to assess performance on a more realistic set of input files.
We thus obtained VCFs from a single-cell RNA-seq study conducted on lung adenocarcinoma patient samples [@Maynard:2020].
The carcinoma VCFs are much smaller, on the order of megabytes rather than gigabytes. 
Alignment was done with STAR and variant calling was performed with GATK HaplotypeCaller. 
The results are shown in *Figure 2* -- `cerebra` clocks in at 34 minutes for the set of 100 VCFs.

[todo: add Figure 2]

One interesting observation is that the first ~10 minutes of the `cerebra` timecourse appear flat, that is, no VCFs are processed.
This can be attributed to the _genome interval tree_ construction phase. 
After the tree is built, files are processed in a near-linear manner. 
Also of note is that `cerebra`'s search operations take advantage of multiprocessing.
Thus `cerebra` should scale better to high-memory machines with more cores, though it has been designed to run on everyday hardware. 

## Conclusions

Our tool satisfies an unmet need in the bioinformatics community. 
Fast and accurate peptide-level summarizing of variants following a sequencing experiment is often crucial to understanding the underlying biology of an experimental system. 
As sequencing costs continue to drop, large-scale variant calling will become accessible to more members of the community, and summary tools like `cerebra` will become increasingly important. 
`cerebra` is fast and accurate and is one of the only tools that fills this niche. 
It offers the advantages of parallel processing and a single, easy-to-interpret output file (CSV or JSON), making downstream analysis accessible to non-bioinformatically inclined members of the community.

## Acknowledgments

Funding for this work was provided by the [Chan Zuckerberg Biohub](https://www.czbiohub.org/). The authors would like
to thank Ashley Maynard, Angela Pisco and Daniel Le for helpful discussions and support.

## Correspondence

Please contact `lincoln.harris@czbiohub.org`

## Code

`cerebra` is written in python. 
Code and detailed installation instructions can be found at https://github.com/czbiohub/cerebra. 
In addition `cerebra` can be found on [PyPi](https://pypi.org/project/cerebra/).

## References

