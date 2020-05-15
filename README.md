cerebra
================================
<a href="https://pypi.org/project/cerebra/"><img alt="PyPI" src="https://badge.fury.io/py/cerebra.svg"></a>

[![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra)
[![Code Coverage](https://codecov.io/gh/czbiohub/cerebra/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/cerebra)

What is _cerebra_?
-------------------------------------

This tool allows you to quickly extract meaningful variant information from a DNA or RNA sequencing experiment. If you're interested in learning what mutations are present in your DNA/RNA samples, variant callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) can be used to generate variant calling format (.vcf) files following a sequencing experiment. However, a single sequencing run can generate on the order of 10^8 unique vcf entries, only a small portion of which contain meaningful biological signal. Thus drawing conclusions from .vcf files remains a substantial challange. _cerebra_ provides a fast and intuitive framework for summarizing vcf entries across samples. It is comprised of four modules that do the following:      

        1) remove germline mutations from samples of interest        
        2) count the total number of mutations in a given sample           
        3) report amino acid level SNPs and indels for each sample             
        4) report the ratio of total to variant reads to each mutation site      
        
_cerebra_ gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

If you're working with tumor data, it might be a good idea to limit the mutational search space to only known cancer variants. Therefore _cerebra_ implements an optional method for restricting to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database.  

NOTE: this framework was developed for, but is certainly not limited to, single-cell RNA sequencing data. 

-   Free software: MIT license


What makes _cerebra_ different from traditional vcf parsers? 
-------------------------------------
Python libraries exist (_ie._ [PyVCF](https://pyvcf.readthedocs.io/en/latest/) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/index.html)) for extracting information from vcf files, and GATK has its own tool for the [task](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php). In fact we integrate vcfpy into our tool. What makes _cerebra_ different is that it reports the RNA transcript and amino acid change associated with each variant. GATK _VariantsToTable_ produces a file that looks like: 
 
    CHROM    POS ID      QUAL    AC
     1        10  .       50      1
     1        20  rs10    99      10

Such a table contains only genomic (ie. DNA-level) coordinates. Often the next question is what specific gene and protein-level mutation is each variant associated with? _cerebra_ queries a reference genome (.fa) and annotation (.gtf) to match each DNA-level variant with its associated gene, probable transcript and probable amino-acid level mutation. _cerebra_ produces a table that looks like the following: 
![alt text](https://raw.githubusercontent.com/lincoln-harris/cerebra/master/cerebra_out_sample.png)

_cerebra_ adheres to HGVS sequence variant [nomenclature](https://varnomen.hgvs.org/) in reporting peptide level variants

Installation
------------

To install the latest version from PyPi you'll first need to install a few system-specific dependencies.     

For OSX:     
```
sudo pip install setuptools
brew update
brew install openssl
brew install zlib
```

For Debian/Ubuntu:     
```
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

Following that, you can install directly from PyPi.        
```pip install cerebra```

If you prefer working with virtual environments you can clone from github and install with `pip`. 
``` 
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
conda create -n cerebra python=3.7
conda activate cerebra
pip install -e . 
````

Usage
-----

_cerebra_ should now be installed as a commandline executable. 
`$ cerebra` should return help information

```
Usage: cerebra  <command>

  high-throughput summarizing of vcf entries following a sequencing
  experiment

Options:
  -h, --help  Show this message and exit.

Commands:
  count-mutations    count total number of mutations in each sample
  find-aa-mutations  report amino-acid level SNPs and indels in each sample
  germline-filter    filter out common SNPs/indels between germline samples...
```


Features
--------
***count-mutations:*** count total number of mutations in each sample         
***find-aa-mutations:*** report amino-acid level SNPs and indels in each sample            
***germline-filter:*** filter out common SNPs/indels between germline samples and samples of interest          


Authors
--------
This work was produced by [Lincoln Harris](https://github.com/lincoln-harris), [Rohan Vanheusden](https://github.com/rvanheusden), [Olga Botvinnik](https://github.com/olgabot) and [Spyros Darmanis](https://spyrosdarmanis.wixsite.com/mylab) of the Chan Zuckerberg Biohub. For questions please contact lincoln.harris@czbiohub.org


Contributing
--------
We welcome any bug reports, feature requests or other contributions. Please submit a well documented report on our [issue tracker](https://github.com/czbiohub/cerebra/issues). For substantial changes please fork this repo and submit a pull request for review. 

Feel free to clone but NOTE this project is still a work in progress. 
