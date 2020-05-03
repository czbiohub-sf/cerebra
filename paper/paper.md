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
