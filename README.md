cerebra
================================

[![image](https://img.shields.io/pypi/v/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://pypi.python.org/pypi/%7B%7B%20cookiecutter.repo_name%20%7D%7D) [![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra)


What is _cerebra_?
-------------------------------------

This tool allows you to extract meaningful variant calls from a single-cell RNA-seq experiment. Mutation callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) generate 1000s of output files (.vcf) for the typical RNA-seq experiment, and consolidating that output into a form from which meaningful biological conclusions can be drawn represents a significant challange. _cerebra_ is comprised of several modules which do the following: 1) generate a cell_x_gene mutation-counts matrix, 2) generate a cell_x_ROI summary table that reports amino acid level mutations for a user-defined list of genes, 3) report read coverage (variant vs total reads) to each ROI.   

-   Free software: MIT license
-   Documentation: <https://>czbiohub.github.io/cerebra

Installation
------------

To install this package, clone this github repository and use pip to install

```
git clone <https://github.com/>czbiohub/cerebra.git 
cd cerebra 

# create a new conda environment
conda create -n cerebra python=3.6
conda activate cerebra

# The '.' means 'install *this*, the folder where I am now'. The 'e' means...?
pip install -e . 
```

Usage
-----

_cerebra_ should now be installed as a commandline executable. 
`$ cerebra` should return help information

```
$ cerebra
Usage: cerebra  <command>

  finds mutants in your scRNA-seq experiment

Options:
  -h, --help  Show this message and exit.

Commands:
  check_coverage_loci           evaluate coverage for each loci for which
                                we...
  check_coverage_whole_gene     evaluate the coverage across every loci
                                with...
  count-mutations
  fusion_search                 searches STAR-fusion output files for a...
  fusions_x_cell                create a cell-wise fusion counts table.
  generate_summary_tables       generate by cell and by sample summary...
  generate_summary_tables_test  generate by cell and by sample summary...
  germline-filter               Given a set of single-cell VCFs and bulk...
  get_aa_mutations              for a specific gene of interest, get the...
  get_mutationalburden          returns the total number of mutations...
  hello                         helloWorld
  s3_import                     import necessary files from s3
```


Features
--------

-   test

