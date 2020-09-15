cerebra
================================
<a href="https://pypi.org/project/cerebra/"><img alt="PyPI" src="https://badge.fury.io/py/cerebra.svg"></a>
[![Docker Build](https://img.shields.io/docker/cloud/build/lincolnharris/cerebra)](https://hub.docker.com/r/lincolnharris/cerebra)    
[![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra) 
[![codecov](https://codecov.io/gh/czbiohub/cerebra/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/cerebra)


What is `cerebra`?
-------------------------------------
This tool allows you to quickly summarize and interpret VCF files.      

If you're interested in learning the full complement of genomic variants present in your DNA/RNA samples, a typical workflow might involve sequencing, followed by variant calling with a tool such as GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php), which generates [variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (VCF) files.       

A VCF file looks like this:

```##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
chr1 631391 . C T 72.28 . AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;QD=25.36;SOR=2.303 GT:AD:DP:GQ:PL 1/1:0,2:2:6:84,6,0
```
Note that only a single VCF record is displayed here.
A sequencing run can generate on the order of 10<sup>8</sup> unique VCF records. 
Only a small portion of these VCF records are likely to be relevant to the researcher.
Thus drawing conclusions from VCF files remains a substantial challenge.


`cerebra` provides a fast and intuitive framework for summarizing VCF records across samples.
It comprises three modules that do the following:      

        1) remove germline variants from samples of interest        
        2) count the total number of variants in a given sample, on a per-gene basis           
        3) report protein variants for each sample                 
        
`cerebra` gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

If you're working with tumor data, it might be a good idea to limit the variant search space to only known, and potentially actionable, cancer variants. 
Therefore `cerebra` implements an optional method for restricting to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database.  

This tool was developed for, but is certainly not limited to, single-cell RNA sequencing data. 
`cerebra` is free software available under the MIT license.


What makes `cerebra` different from traditional VCF parsers? 
-------------------------------------
Python libraries exist (_i.e._ [PyVCF](https://pyvcf.readthedocs.io/en/latest/) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/index.html)) for extracting information from VCF files.
Another is [GATK VariantsToTable](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php), which produces a file that looks like this:
 
    CHROM    POS ID      QUAL    AC
     1        10  .       50      1
     1        20  rs10    99      10

This table contains only genomic (_i.e._ DNA-level) coordinates. 
Often the next questions are: what are the genes associated with these variants, and what are the peptide-level effects of these variants?
`cerebra` queries a reference genome (.fa) and annotation (.gtf) to match each DNA-level variant with its associated gene, and its probable protein variant.
`cerebra` produces the following outfile: 

```
$ python
> import json
> f = open(/path/to/cerebra/output.json)
> data = json.load(f)
> print(json.dumps(data, indent=4))

{
    "CCN1": {
        "A16_B000563": [],
        "A1_B001546": [],
        "A1_B002531": [
            "ENSP00000398736.2:p.(Glu189=)"
        ],
        "A1_B002570": [],
        "A2_B002558": [],
        "A3_B000561": [
            "ENSP00000398736.2:p.(Arg209Trp)",
            "ENSP00000398736.2:p.(Ile90Val)"
        ],
        "A3_B000568": [],
        "A3_B001544": [
            "ENSP00000398736.2:p.(Ala82Thr)"
        ],
        "A3_B002090": [],
        "A3_B002531": [
            "ENSP00000398736.2:p.(Pro217Ser)"
        ]
    }
}
```

Here _CCN1_ is a gene name while _A16_B000563_, _A1_B001546_, _A1_B002531,_... are RNA-seq sample IDs.
`cerebra` reports variants to every gene in the genome, for every sample in a given experiment. 
The _ENSP*_ numbers refer to [Ensembl](https://www.ensembl.org/index.html) translation IDs -- unique identifiers that correspond to exactly one polypeptide in the Ensembl database. 
The strings enclosed in parentheses refer to the amino-acid level variants reported in that particular sample. 
For example the string `Arg209Trp` indicates that position 209 of this particular polypeptide should contain an _Arg_, but the experimental sample instead contains a _Trp_. 
`cerebra` adheres to HGVS sequence variant [nomenclature](https://varnomen.hgvs.org/) in reporting amino-acid variants.


Features
--------
### `germline-filter`

Variant calling is often applied to the study of cancer. 
If the research project is centered around a “tumor vs. normal” question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between tumor and normal samples, and thus excludes variants unlikely to be pathogenic for the cancer under study.
The user provides a very simple metadata file (see [USAGE.md](https://github.com/czbiohub/cerebra/blob/master/docs/USAGE.md)) that indicates which tumor samples correspond to which normal samples.
The output of germline-filter is a set of trimmed-down VCF files, which will be used for the next two steps.
If you do not have access to “normal” samples then proceed directly to `count-variants` or `find-peptide-variants`.

### `count-variants`
This module reports the raw variant counts for every gene across every sample.
The output is a CSV file that contains counts for each sample versus every gene in the genome. 
If working with cancer variants, the user has the option of limiting the search space to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database. 

### `find-peptide-variants`
This module reports the protein variations associated with each genomic variant.
VCF records are converted to peptide-level variants, and then [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, 
in accordance to the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature. 
The output is a hierarchically ordered text file (CSV or JSON) that reports the Ensemble protein ID and the gene associated with each variant, for each experimental sample. 
The user again has the option of limiting the search space to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database. 
This module also has an option to report the number of variant vs. wildtype reads found at each locus. 


Dependencies
------------
`cerebra` depends on some (fairly standard) packages and libraries. 
Before installing, it might be a good idea to make sure all of the requisite packages are installed on your system (_note:_ if installing with Docker you can skip this step). 

__MacOS Dependencies:__
```
sudo pip install setuptools
brew update
brew install openssl
brew install zlib
```

__Linux Dependencies:__
```
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

As of present `cerebra` is not installable on Windows. 
`cerebra` depends on the [`pysam`](https://pysam.readthedocs.io/en/latest/index.html) library (or rather, `pysam` is a dependency-of-a-dependency) and currently this library is only available on Unix-like systems. 
Windows solutions like [WSL](https://docs.microsoft.com/en-us/windows/wsl/) exist for overcoming precisely this challenge, however, `cerebra` has not been tested on WSL or any other Unix-like subsystem for Windows.    


Installation (for users)
------------
There are four different methods available to install `cerebra`.
Choose one of the following:

__With [Docker](https://hub.docker.com/r/lincolnharris/cerebra) (recommended)__          
```
docker pull lincolnharris/cerebra
docker run -it lincolnharris/cerebra
```      
_warning: this image will take up ~1Gb of space._               

__With the python standard library [`venv`](https://docs.python.org/3/library/venv.html) module__
```
python3 -m venv cerebra
source cerebra/bin/activate
pip3 install cerebra 
```

__With [conda](https://docs.conda.io/en/latest/)__
``` 
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install cerebra
```

__From [PyPi](https://pypi.org/project/cerebra/) (system-wide installation, NOT RECOMMENDED)__    
For novice users, it's generally a better idea to install packages within virtual environments. 
However, `cerebra` can be installed system-wide, if you so choose. 
```
pip3 install cerebra

# OR, if you dont have root privileges
pip3 install --user cerebra
```


Installation (for developers)
------------
Here's how to set up cerebra for local development. 
After installing the requisite [dependencies](https://github.com/czbiohub/cerebra/blob/master/README.md#dependencies):

1.  Fork the `cerebra` repo on GitHub: https://github.com/czbiohub/cerebra
2.  Clone your fork locally:

        $ git clone https://github.com/your-name/cerebra.git

3.  Install your local copy into a virtualenv. Using the standard library [`venv`](https://docs.python.org/3/library/venv.html) module: 

        $ cd cerebra
        $ python3 -m venv cerebra-dev
        $ source cerebra-dev/bin/activate
        $ pip3 install -r requirements.txt -r test_requirements.txt -e .

4.  Create a branch for local development:

        $ git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5.  When you're done making changes, check that your changes pass `flake8` and `pytest`:

        $ make test
        $ make coverage
        $ make lint

6.  Commit your changes and push your branch to GitHub:

        $ git add .
        $ git commit -m "Your detailed description of your changes."
        $ git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.
See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for more. 


(Quickstart) Usage
-----
`$ cerebra` should return usage information

```
Usage: cerebra  <command>

  a tool for fast and accurate summarizing of variant calling format (VCF)
  files

Options:
  -h, --help  Show this message and exit.

Commands:
  germline-filter    filter out common SNPs/indels between tumor and normal samples
  count-variants    count total number of variants in each sample, and report on a per-gene basis
  find-peptide-variants  report peptide-level SNPs and indels in each sample, and associated coverage
```

Note that the `-h` command will display usage information for each of the three commands. 

An example workflow might look like this:   

**Step 1:**     
```
cerebra germline-filter --processes 2 --normal_path /path/to/normal/vcfs --tumor_path /path/to/tumor/vcfs --metadata /path/to/metadata/file --outdir /path/to/filtered/vcfs
```

**Step 2:**     
```
cerebra count-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --refgenome /path/to/genome/annotation --outfile /path/to/output/file /path/to/filtered/vcfs/*
```

**Step 3:**          
```
cerebra find-peptide-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --annotation /path/to/genome/annotation --genomefa /path/to/genome/fasta --report_coverage 1 --output /path/to/output/file /path/to/filtered/vcfs/*
```

For advanced usage information, see [USAGE.md](https://github.com/czbiohub/cerebra/blob/master/docs/USAGE.md). 


Authors
--------
This work was produced by [Lincoln Harris](https://github.com/lincoln-harris), [Rohan Vanheusden](https://github.com/rvanheusden), [Olga Botvinnik](https://github.com/olgabot) and [Spyros Darmanis](https://spyrosdarmanis.wixsite.com/mylab) of the Chan Zuckerberg Biohub. 
For questions please contact ljharris018@gmail.com


Contributing
--------
We welcome any bug reports, feature requests or other contributions. 
Please submit a well documented report on our [issue tracker](https://github.com/czbiohub/cerebra/issues). 
For substantial changes please fork this repo and submit a pull request for review. 

See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for additional details. 

You can find official releases [here](https://github.com/czbiohub/cerebra/releases). 
