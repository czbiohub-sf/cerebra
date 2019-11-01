cerebra
================================

[![image](https://img.shields.io/pypi/v/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://pypi.python.org/pypi/%7B%7B%20cookiecutter.repo_name%20%7D%7D)


What is _cerebra_?
-------------------------------------

This tool allows you to extract meaningful variant calls from a single-cell RNA-seq experiment. Mutation callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) generate 1000s of output files (.vcf) for a typical single-cell RNA-seq experiment, and consolidating that output into a form from which meaningful biological conclusions can be drawn represents a significant challange. _cerebra_ is comprised of several modules which do the following: **1)** generate a cell_x_gene mutation-counts matrix, **2)** generate a cell_x_ROI summary table that reports amino acid level mutations for a user-defined list of genes, **3)** report read coverage (variant vs total reads) to each ROI. _cerebra_ gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

-   Free software: MIT license
-   Documentation: <https://czbiohub.github.io/cerebra>

Installation
------------

To install this package, clone this github repository and use pip to install

```
git clone <https://github.com/>czbiohub/cerebra.git 
cd cerebra 

# create a new conda environment
conda create -n cerebra python=3.6
conda activate cerebra

# The '.' means install *this*, the folder where I am now.
pip install -e . 
```

Usage
-----

_cerebra_ should now be installed as a commandline executable. 
`$ cerebra` should return help information

```
Usage: cerebra  <command>

  finds mutants in your scRNA-seq experiment

Options:
  -h, --help  Show this message and exit.

Commands:
  count-mutations    count total number of mutations in each sample
  find-aa-mutations  report amino-acid level SNPs and indels in each sample
  germline-filter    filter out common SNPs/indels between germline samples...
  get-coverage       report coverage to each SNP location contained within
                     a...
```


Features
--------
***count-mutations:*** count total number of mutations in each sample         
***find-aa-mutations:*** report amino-acid level SNPs and indels in each sample            
***germline-filter:*** filter out common SNPs/indels between germline samples and samples of interest          
***get-coverage:*** report coverage to each SNP location contained within a set of genes  


Authors
--------
This work was produced by Lincoln Harris, Rohan Vanheusden and Spyros Darmanis of the Chan Zuckerberg Biohub


Correspondence
--------
For questions please contact lincoln.harris@czbiohub.org

