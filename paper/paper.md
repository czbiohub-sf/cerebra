---
title: '`cerebra`: A tool for fast and accurate summarizing of variant calling format (VCF) files'
tags:
    - python
    - genomics
    - variant calling
    - VCF
    - single-cell
    - cancer
authors:
 - name: Lincoln Harris
   orcid: 0000-0003-2544-4225
   affiliation: 1
 - name: Rohan Vanheusden
   affiliation: 1 
 - name: Olga Botvinnik
   orcid: 0000-0003-4412-7970
   affiliation: 1 
 - name: Spyros Darmanis
   orcid: 0000-0003-4002-8158
   affiliation: 1 
affiliations: 
- name: Chan Zuckerberg Biohub, San Francisco, CA
  index: 1
date: 20, August, 2020 
bibliography: paper.bib
---

## Motivation

A single "typo" in the genome can have profound consequences on an organism's biology.
Identifying the protein changes that accompany these genomic typos is a fundamental challenge in bioinformatics. 
Tools exist for identifying genomic variants and predicting their associated peptide-level changes, however, wrangling genomic variant calls and peptide-level predictions across thousands of samples remains an substantial challenge. 
`cerebra` addresses this need by offering a fast and accurate framework for summarizing genomic variant calls and peptide-level predictions across many samples. 

To find variants in the genome, researchers often begin with a [DNA-sequencing](https://en.wikipedia.org/wiki/DNA_sequencing) (DNA-seq) or [RNA-sequencing](https://en.wikipedia.org/wiki/RNA-Seq) (RNA-seq) experiment on their samples of interest.
Sequencing is followed by alignment of reads to the reference genome with tools like [STAR](https://github.com/alexdobin/STAR) or [BWA](http://bio-bwa.sourceforge.net/), followed by variant calling with tools like [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) 
or [freebayes](https://github.com/ekg/freebayes) [@star; @bwa; @haplocaller; @freebayes]. 
Variant callers produce tab delimited text files in the [variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (VCF) for each processed sample.
VCF files encode: i) the genomic position, ii) reference vs. observed DNA sequence, and iii) quality associated with each observed variant. 
Shown below are the first 4 lines of a sample VCF file. 
Note that only a single record is displayed, and that the record line has been artificially wrapped.

```
##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM  POS ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
chr1	631391	.	C	T	72.28	.	AC=2;AF=1.00;AN=2;DP=2;
        ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;
        QD=25.36;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:84,6,0
```

Current methods for variant calling are incredibly powerful and robust, however, a single sequencing run can generate as many as $10^{8}$ unique VCF records.
Only a small portion of these VCF records are likely to be relevant to the researcher.
In addition, variant callers report only the genomic location and not the likely effect the variant has on the translated protein sequence.
To address the unmet need for high-throughput VCF summary tools, we introduce `cerebra`, a python package that provides fast and accurate peptide-level summarizing of VCF files.


## Functionality

`cerebra` comprises three modules: i) `germline-filter` removes variants that are common between control/germline samples 
and samples of interest, ii) `count-variants` reports total number of variants in each sample, and iii) `find-peptide-variants` reports likely peptide-level variants in each sample. 
Here we use _variant_ to refer to single nucleotide polymorphisms (SNPs) and short insertions/deletions. 
`cerebra` is not capable of reporting larger structural variants such as copy number variations and chromosomal rearrangements.

A data structure crucial to `cerebra` is the *genome interval tree*, which matches RNA transcripts
and peptides to each feature in the genome (\autoref{workflow}). 
[Interval trees](https://en.wikipedia.org/wiki/Interval_tree) are self-balancing binary search trees that store numeric intervals and can quickly retrieve every such interval that overlaps a given query interval (_[see also](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_). 
Given _n_ nodes, interval trees have theoretical average-case O(log*n*) and worst-case O(*n*) time complexity for search operations, making them tractable for genome-scale operations [@Cormen:2009, _[see also](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_].
Tree construction proceeds at O(*n*log*n*) time complexity, making construction rather than search the bottleneck for most VCF sets [@Alekseyenko:2007]. 
The _genome interval tree_ is constructed with a reference genome sequence ([FASTA format](https://en.wikipedia.org/wiki/FASTA_format), often with a `.fa` extension), and a genome annotation 
([gene transfer format, GTF](https://www.gencodegenes.org/pages/data_format.html), `.gtf` extension).
We rely on the [ncls](https://github.com/biocore-ntnu/ncls) python library for fast interval tree construction and lookup operations [@Alekseyenko:2007].

In order to analyze multiple VCF records at once, we use [multiprocessing](https://en.wikipedia.org/wiki/Multiprocessing) with the Python [pathos](https://pypi.org/project/pathos/) library module [@pathos].
We extract relevant information -- including genomic interval, observed base, and read coverage -- from each variant record. 
In the `germline-filter` module variants are compared to one another and filtered out if found to be identical.
In `count-variants` variants are simply matched to whichever gene they came from. 
In `find-peptide-variants` variants are queried against our _genome interval tree_ -- if a matching interval is found we convert the DNA-level variant to a peptide-level variant. 
Finally, peptide-level variants across all VCF files are reported in tabular format. 

![Workflow describing the `find-peptide-variants` module.
We construct a genome interval tree from a genome annotation (.gtf) and a reference genome sequence (.fa), then process VCF files in parallel to create a single tabular output file (CSV or JSON).\label{workflow}](fig1.jpg)

## `germline-filter`

If the research project is centered around a "tumor/pathogenic vs control" question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between the control and the experimental tissue so as to not bias the results by including non-pathogenic variants. 
The user provides a very simple metadata file (see [README.md](https://github.com/czbiohub/cerebra/blob/master/README.md)) that indicates which experimental samples correspond to which control samples.
Using the [vcfpy](https://pypi.org/project/vcfpy/) library we quickly identify shared variants across control/experimental matched VCF files, then write new VCFs that contain only the unique variants [@vcfpy].
These steps are performed by a [subprocess pool](https://pypi.org/project/pathos/) so that we can process multiple discreet chunks of input at the same time. 

The output of `germline-filter` is a set of trimmed-down VCF files, which will be used for the next two steps. 
If you do not have access to "control" tissue then proceed directly to `count-variants` or `find-peptide-variants`. 

## `count-variants`

The `count-variants` module reports the raw variant counts for every gene across every sample.
We first create a _genome interval tree_ from the reference GTF, then read in a VCF file and convert it to a [vcfpy](https://pypi.org/project/vcfpy/) object, then processes VCF records in [parallel](https://en.wikipedia.org/wiki/Multiprocessing). 
Each variant is matched to its corresponding gene, and gene-wise counts are stored in [shared memory](https://en.wikipedia.org/wiki/Shared_memory).

If working with cancer samples, the user has the option to limit the reported variants to those also found in Wellcome Sanger Institute's [COSMIC](https://cancer.sanger.ac.uk/cosmic) database [@cosmic]. 
While certainly not exhaustive, this database contains an extensive list of known human variants. 
This option is designed to limit the search space to known and potentially actionable targets. 

`count-variants` produces two output files, one containing raw variant counts and one containing COSMIC filtered variant counts for every gene in the genome. 

## `find-peptide-variants`

The `find-peptide-variants` module reports the peptide-level consequence of genomic variants.
First we load the reference GTF, then construct an index (.fai) of the genome fasta file with [pyfaidx](https://pypi.org/project/pyfaidx/) to enable fast random memory access [@pyfaidx].
We then create a _genome interval tree_ that will be used to quickly match genomic coordinates from VCF records to peptide-level variants. 
The user again has the option to limit the search space to variants found in the COSMIC database. 
VCF records are read in simultaneously; individual records are converted to _GenomePosition_ objects to keep track of their genomic intervals and observed DNA bases.
_GenomePositions_ are then queried against the _genome interval tree_. 
If an overlapping interval is found we retrieve the peptide-level variant from this node of the _genome interval tree_. 
Peptide-level variants are converted to [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, in accordance with the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature [@ensembl; @hgvs].
The output is a hierarchically ordered text file (CSV or JSON) that reports the the ENSEMBL protein ID and the gene associated with each variant, for each experimental sample.    

Variant callers are known to produce a great deal of false positives, especially when applied to single-cell RNA-seq data [@Enge:2017].
To address this concern we have included the `--report_coverage` option. 
If indicated this option will report counts for both variant and wildtype reads at all variant loci. 
We reasoned that variants with a high degree of read support are less likely to be false positives.
This option is designed to give the user more confidence in individual variant calls.        

We should emphasize that `find-peptide-variants` does not *definitively* report peptide-level variants but rather the *likely* set of peptide variants. 
Definitively reporting protein variants from RNA-seq requires knowledge of alternate splicing -- this represents an open problem in the field [@Huang:2017]. 
For example, if a read picks up a variant in exon 2 of a given gene, we can report each of the potential spliceforms of that gene that contain exon 2, but we **cannot** infer which of those particular spliceforms are actually present in our sample (see \autoref{splice}). 
For the example shown in \autoref{splice} we would translate and report _t1_ and _t3_ as both of these contain exon 2. 
It is possible the sample does not actually express both of these spliceforms, however, determining the spliceform landscape of a sample from RNA-seq is outside the scope of this project. 

![For a given mutational event, `cerebra` reports ALL potentially affected spliceforms.\label{splice}](fig2.jpg)

To assess performance of `find-peptide-variants` we obtained VCFs from a single-cell RNA-seq study conducted on lung adenocarcinoma patient samples [@Maynard:2019]. 
These VCFs were produced with STAR (alignment) and GATK HaplotypeCaller (variant calling), and are on the order of megabytes, typical of a single-cell RNA-seq experiment. 
`cerebra` was run on standard hardware (MacBook Pro, 2.5GHz quad-core processor, 16 GB RAM).
As show in \autoref{runtime} `cerebra` processed a set of 100 VCF files in approximately 34 minutes. 

![`cerebra` processes 100 VCF files (~400 Mb in total) in ~34 minutes.\label{runtime}](fig3.jpg)

The first 10 or so minutes of `cerebra find-peptide-variants` do not involve any VCF processing, instead, this time is attributed to the _genome interval tree_ construction phase.
After the tree is built, files are processed in a near-linear manner. 
Also of note is that `cerebra`'s search operations take advantage of multiprocessing.
`cerebra` should scale better to high-memory machines than single-threaded tools, though it has been designed to run on standard hardware.

## Conclusions

RNA/DNA sequencing paired with fast and accurate summarizing of variants is often crucial to understanding the biology of an experimental system. 
We present a tool that can be used to quickly summarize the variant calls contained within a large set of VCF files.
As sequencing costs continue to drop, large-scale variant calling will become accessible to more members of the community, and summary tools like `cerebra` will become increasingly important. 
Our tool offers the advantages of parallel processing and a single, easy-to-interpret output file (CSV or JSON).

`cerebra` is already enabling research, see [@Maynard:2019], a study that examines the tumor microenvironment of late-stage drug-resistant carcinomas with single-cell RNA-sequencing. 
Understanding the mutational landscape of individual tumors was essential to this study, and given the sheer volume of VCF records, would not have been possible without `cerebra`. 
We hope that `cerebra` can provide an easy-to-use framework for future studies in the same vein. 

## Acknowledgments

Funding for this work was provided by the [Chan Zuckerberg Biohub](https://www.czbiohub.org/).
The authors would like to thank Ashley Maynard, Angela Pisco and Daniel Le for helpful discussions and support.

## Correspondence

Please contact `ljharris018@gmail.com`

## Code

`cerebra` is written in Python 3. 
Code and detailed installation instructions can be found at https://github.com/czbiohub/cerebra. 
In addition `cerebra` can be found on [PyPi](https://pypi.org/project/cerebra/).

## References

