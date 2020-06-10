`cerebra`
================================
<a href="https://pypi.org/project/cerebra/"><img alt="PyPI" src="https://badge.fury.io/py/cerebra.svg"></a>

[![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra)
[![Code Coverage](https://codecov.io/gh/czbiohub/cerebra/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/cerebra)

What is `cerebra`?
-------------------------------------

This tool allows you to quickly extract meaningful variant information from a DNA or RNA sequencing experiment. If you're interested in learning what variants are present in your DNA/RNA samples, variant callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) can be used to generate variant calling format (VCF) files following a sequencing experiment. However, a single sequencing run can generate on the order of 10^8 unique VCF entries, only a small portion of which contain meaningful biological signal. Thus drawing conclusions from VCF files remains a substantial challange. `cerebra` provides a fast and intuitive framework for summarizing VCF entries across samples. It is comprised of three modules that do the following:      

        1) remove germline variants from samples of interest        
        2) count the total number of variants in a given sample           
        3) report peptide-level variants for each sample                 
        
`cerebra` gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

If you're working with tumor data, it might be a good idea to limit the variant search space to only known cancer variants. Therefore `cerebra` implements an optional method for restricting to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database.  

NOTE: this framework was developed for, but is certainly not limited to, single-cell RNA sequencing data. 

-   Free software: MIT license


What makes `cerebra` different from traditional VCF parsers? 
-------------------------------------
Python libraries exist (_ie._ [PyVCF](https://pyvcf.readthedocs.io/en/latest/) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/index.html)) for extracting information from VCF files, and GATK has its own tool for the [task](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php). In fact we integrate vcfpy into our tool. What makes `cerebra` different is that it reports the RNA transcript and amino acid change associated with each variant. GATK _VariantsToTable_ produces a file that looks like: 
 
    CHROM    POS ID      QUAL    AC
     1        10  .       50      1
     1        20  rs10    99      10

Such a table contains only genomic (ie. DNA-level) coordinates. Often the next question is what specific gene and peptide-level variants is each variant associated with? `cerebra` queries a reference genome (.fa) and annotation (.gtf) to match each DNA-level variant with its associated gene, probable transcript and probable peptide-level level variants. `cerebra` produces a table that looks like the following: 
![alt text](https://raw.githubusercontent.com/lincoln-harris/cerebra/master/cerebra_out_sample.png)

`cerebra` adheres to HGVS sequence variant [nomenclature](https://varnomen.hgvs.org/) in reporting peptide level variants

Features
--------
### `germline-filter`

If the research project is centered around a "tumor/pathogenic vs control" question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between the control and the experimental tissue so as to not bias the results by including non-pathogenic variants. 
The user provides a very simple metadata file that indicates which experimental samples correspond to which control samples. For example:

```
experimental_sample_id,germline_sample_id
sample1,gl_sample1
sample2,gl_sample1
sample3,gl_sample2
sample4,gl_sample2
sample5,gl_sample2
```

There is also the option to limit the reported variants to those found in NCBI's [dbSNP](https://www.ncbi.nlm.nih.gov/books/NBK21088/) and the Wellcome Sanger Institute's [COSMIC](https://cancer.sanger.ac.uk/cosmic) databases. 
This option is designed to give the user a higher degree of confidence in the pathogenic nature of each variant -- if independent experiments have reported a given variant in human tissue, there is a higher likilihood that it is pathogenic. 
The output of `germline-filter` is a set of trimmed-down VCF files. 

If you have access to "control" tissue and your experimental question is concerned with differences between tumor/pathogenic tissue and control tissue, then `germline-filter` is the right place to start.
`germline-filter` will produce a new set of VCFs, which you'll use for the next two steps.
If you do not have access to "control" tissue, then proceed directly to `count-mutations` or `find-aa-mutations`.

### `count-mutations`
The `count-mutations` module reports the raw variant counts for every gene across every sample.
The output is a CSV file that contains counts for each sample versus every gene in the genome. 

### `find-aa-mutations`
The `find-aa-mutations` module reports the peptide-level consequence of variants in the genome.
If working  with cancer samples, the user has the option to filter out all variants that are not found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database and are therefore unlikely to be pathogenic.
VCF records are converted to peptide-level variants, and then [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, 
in acordance to the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature. 
The output is a heirarchically ordered text file (CSV or JSON) that reports the the Ensemble protein ID and the gene associated with each variant, for each experimental sample. 

We should stress that `find-aa-mutations` does not *definitively* report peptide-level variants but rather the *likely*
set of peptide variants. 
Definitively reporting protein variants requires knowledge of alternate splicing -- this represents an open problem in scRNA-seq.
For example, if a read picks up a variant in exon 2 of geneA, we can report each of the potential spliceforms of geneA that contain exon 2, but we **cannot** infer which of those particular spliceforms are actually present in our sample. 
Thus we report all possible spliceforms; determining the spliceform landscape of an individual cell from scRNA-seq is outside the scope of this project. 


Note that all three modules take advantage of multiprocessing.
Thus `cerebra` should scale better to high-memory machines with more cores, though it has been designed to run on everyday hardware. 

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

`cerebra` should now be installed as a commandline executable. 
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

Note that the `-h` command will display usage information for each of the three commands. 

An example workflow might look like this:   

**Step 1:**     
`cerebra germline-filter --processes 2 --germline /path/to/control/vcfs --cells /path/to/experimental/vcfs --metadata /path/to/metadata/file --outdir /path/to/filtered/vcfs`   

**Step 2:**     
`cerebra count-mutations --processes 2 --cosmicdb /optional/path/to/cosmic/database --refgenome /path/to/genome/annotation --outfile /path/to/output/file /path/to/filtered/vcfs/*` 

**Step 3:**          
`cerebra find-aa-mutations --processes 2 --cosmicdb /optional/path/to/cosmic/database --annotation /path/to/genome/annotation --genomefa /path/to/genome/fasta --report_coverage 1 --output /path/to/output/file /path/to/filtered/vcfs/*`


Authors
--------
This work was produced by [Lincoln Harris](https://github.com/lincoln-harris), [Rohan Vanheusden](https://github.com/rvanheusden), [Olga Botvinnik](https://github.com/olgabot) and [Spyros Darmanis](https://spyrosdarmanis.wixsite.com/mylab) of the Chan Zuckerberg Biohub. For questions please contact lincoln.harris@czbiohub.org


Contributing
--------
We welcome any bug reports, feature requests or other contributions. Please submit a well documented report on our [issue tracker](https://github.com/czbiohub/cerebra/issues). For substantial changes please fork this repo and submit a pull request for review. 

Feel free to clone but NOTE this project is still a work in progress. 
You can find official releases on [PyPi](https://pypi.org/project/cerebra/)
