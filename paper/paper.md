---
title: 'cerebra: A tool for fast and accurate summarizing of variant calling format (.vcf) files'
authors:
- affiliation: 1
  name: Lincoln Harris
  orcid: 
- affiliation: 1
  name: Rohan Vanheusden
  orcid: 
- affiliation: 1
  name: Olga Borisovna Botvinnik
  orcid: 0000-0003-4412-7970
- affiliation: 1
  name: Spyros Darmanis
  orcid: 

date: "May 2020"
output:
  html_document:
    df_print: paged
  pdf_document: default
  word_document: default
bibliography: paper.bib
tags:
- python
- genomics
- variant calling
- vcf
- single-cell
affiliations:
- index: 1
  name: Chan Zuckerberg Biohub, San Francisco, CA
---


# Motivation
Researchers are often interested in identifying the DNA mutations, or variants, present in a set of samples.
Typically a DNA or RNA sequencing experiment (DNA-seq, RNA-seq) is performed, followed by alignment to the 
reference genome with tools like [STAR](https://github.com/alexdobin/STAR) and 
[BWA](http://bio-bwa.sourceforge.net/), followed by variant calling with tools like [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) 
and [freeBayes](https://github.com/ekg/freebayes). Variant callers produce tab
delimited text files ([variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf), *.vcf*)
for each processed sample, which encode the genomic position, reference vs. observed DNA sequence, and quality
associated with each observed variant. Current methods for variant calling are incredibly powerful and robust, 
however, a single sequencing run can generate on the order of 10^8 unique vcf entries, only a 
small portion of which are of relevance to the researcher. In addition variant callers only report DNA level 
mutations and not the functional consequences of each mutation, *ie.* peptide-level variants. We introduce
*cerebra*, a python package that provides fast and accurate peptide-level summarizing of vcf files.

# Functionality
*cerebra* comprises three modules: **germline-filter** removes variants that are common between germline samples and samples of interest, **count-mutations** reports total number of variants in each sample, and **find-aa-mutations** reports likely peptide-level variants in each sample. Here *variants* refers to single nucleotide polymorphisms (SNPs) and short insertions and deletions. *cerebra* is not capable of reporting larger structural variants such as copy number variations and chromosomal rearrangements. *cerebra* utilizes a genome sequence (.fa) database, a genome annotation (.gtf) and a reference transcriptome to construct a *genome interval tree*, a data structure that matches RNA transcripts and peptides to each feature in the genome (*Figure 1*). Interval trees are self-balancing binary search trees that store numeric intervals and can quickly find every such interval that overlaps a given query interval. They have theoretical average-case O(log*n*) and worst-case O(*n*) time complexity for search operations, making them tractable for genome-scale operations. Tree construction proceeds at O(*n*log*n*) time complexity, making construction rather than search the bottleneck for small vcf sets. [todo: say something about memory footprint]

Recently, a beta (*ie.* development) version of a similar tool was released by the Genome Analysis Toolkit ([GATK](https://software.broadinstitute.org/gatk/)) project of the Broad Institute. We set out to benchmark performance of *cerebra* against this tool, called GATK [Funcotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php). For initial comparison we used a high-quality reference vcf set from the [Genome in a Bottle Consortium](https://www.nist.gov/programs-projects/genome-bottle) (GiaB). This is a set of seven large (~2 GB), reference-grade vcf files from two Mendelian Trios and one individual; variant calls have been verified by independent sequencing platforms, alignment tools and variant callers. *cerebra* and Funcotator were each run on the GiaB vcf set on a laptop (2.5GHz quad-core processor, 16 GB RAM). A comparison of run-times can be found in *Figure 2*; *cerebra* processes the seven files in 44 minutes compared to Funcotator's 195 minutes.

The GiaB vcfs are quite large and are perhaps not comparable to those generated in a typical sequencing experiment. Thus we compared performance on an independent set of vcf files. This vcf set comes from a single-cell RNA-seq experiment \cite{maynard2019lung} conducted on lung adenocarcinoma patient samples. Alignment was done with STAR and variant calling performed with GATK HaplotypeCaller. The carcinoma vcfs are much smaller, on the order of megabytes rather than gigabytes. The results are shown in *Figure 2* -- again, *cerebra* vastly outperforms Funcotator, clocking in at 34 and 139 minutes respectively for the set of 100 vcfs.

One interesting observation is that the first ~10 minutes of the *cerebra* timecourse appear flat, that is, no vcfs are processed. This can be attributed to the genome interval tree construction phase. After the tree is built, files are processed in a near-linear manner. Also of note is that *cerebra*'s search operations take advantage of multiprocessing while Funcotator is limited to a single process. Thus *cerebra* should scale better to high-memory machines with more cores, though it has been designed to run on everyday hardware. 

In terms of accuracy, we report 91\% concordance between the variant call summaries made by *cerebra* and Funcotator. This is to say that on average 91% of the translation IDs reported by *cerebra* are also reported by Funcotator, for the GiaB reference vcf set (*Figure 3*). [todo: more here] 

Another key difference is that Funcotator reports transcript summaries in a new vcf file for each sample, whereas *cerebra* aggregates transcript summaries across all samples into a single tab-delimited text file. Funcotator thus goes from many vcfs to many vcfs, while *cerebra* goes from many vcfs to a single table. For those less bioinformatically inclined it is far easier to draw conclusions from a single table than many difficult to parse vcf files.

[todo: figure 4]
[todo: table 1]

# Conclusion
We believe our tool satisfies an unmet need in the bioinformatics community. Fast and accurate peptide-level summarizing of variants in a sequencing experiment is often crucial to understanding the underlying biology of an experimental system. As sequencing costs continue to drop, large-scale variant calling will become more accessible to all members of the community, and summary tools like *cerebra* will become more important. #cerebra# is significantly faster and comparably accurate to existing tools that fill this niche. It offers the additional advantage of a single easy-to-interpret output file, making downstream variant calling analysis accessible to non-bioinformatically inclined members of the community.

# Acknowledgments
Funding for this work was provided by the [Chan Zuckerberg Biohub](https://www.czbiohub.org/). The authors would like to thank Ashley Maynard, Angela Pisco and Daniel Le for helpful discussions and support.

# Correspondence
Please contact `lincoln.harris@czbiohub.org`

# Code
*cerebra* is written in python. Code and detailed installation instructions can be found at https://github.com/czbiohub/cerebra. In addition *cerebra* can be found on [PyPi](https://pypi.org/project/cerebra/)

# References
todo: add

# Workflow

![checkout](workflow.jpg)

