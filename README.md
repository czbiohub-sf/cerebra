cerebra
================================

[![image](https://img.shields.io/pypi/v/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://pypi.python.org/pypi/%7B%7B%20cookiecutter.repo_name%20%7D%7D)


What is _cerebra_?
-------------------------------------

This tool allows you to quickly extract meaningful variant information from a DNA or RNA sequencing experiment. If you're interested in learning what mutations are present in your DNA/RNA samples, variant callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) can be used to generate variant calling format (.vcf) files following a sequencing experiment. However, a single sequencing run can generate on the order of 10^8 unique vcf entries, only a small portion of which contain meaningful biological signal. Thus drawing conclusions from .vcf files remains a substantial challange. _cerebra_ provides a fast and intuitive framework for summarizing vcf entries across samples. It is comprised of four modules that do the following: **1)** remove germline mutations from samples of interest, **2)** count the total number of mutations in a given sample, **3)** report amino acid level SNPs and indels for each sample, and **4)** report the ratio of total to variant reads to each mutation site. _cerebra_ gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

NOTE: this framework was developed for, but is certainly not limited to, single-cell RNA sequencing data. 

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

