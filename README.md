cerebra
================================
<a href="https://pypi.org/project/cerebra/"><img alt="PyPI" src="https://badge.fury.io/py/cerebra.svg"></a>
[![Docker Build](https://img.shields.io/docker/build/lincolnharris/cerebra)](https://hub.docker.com/r/lincolnharris/cerebra)    
[![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra)
[![Code Coverage](https://codecov.io/gh/czbiohub/cerebra/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/cerebra)        

What is `cerebra`?
-------------------------------------
This tool allows you to quickly extract meaningful variant information from a DNA or RNA sequencing experiment. If you're interested in learning what variants are present in your DNA/RNA samples, variant callers like GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) can be used to generate variant calling format (VCF) files following a sequencing experiment. A VCF file looks like this:

```##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
chr1 631391 . C T 72.28 . AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;QD=25.36;SOR=2.303 GT:AD:DP:GQ:PL 1/1:0,2:2:6:84,6,0
```
Note that only a single VCF record is displayed here.
A sequencing run can generate on the order of 10^8 unique VCF records, only a small portion of which contain meaningful biological signal. Thus drawing conclusions from VCF files remains a substantial challange. `cerebra` provides a fast and intuitive framework for summarizing VCF records across samples. It is comprised of three modules that do the following:      

        1) remove germline variants from samples of interest        
        2) count the total number of variants in a given sample, on a per-gene basis           
        3) report peptide-level variants for each sample                 
        
`cerebra` gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

If you're working with tumor data, it might be a good idea to limit the variant search space to only known cancer variants. Therefore `cerebra` implements an optional method for restricting to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database.  

This tool was developed for, but is certainly not limited to, single-cell RNA sequencing data. 

-   Free software: MIT license


What makes `cerebra` different from traditional VCF parsers? 
-------------------------------------
Python libraries exist (_i.e._ [PyVCF](https://pyvcf.readthedocs.io/en/latest/) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/index.html)) for extracting information from VCF files, and GATK has its own tool for the [task](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php). In fact we integrate vcfpy into our tool. What makes `cerebra` different is that it reports the RNA transcript and amino acid change associated with each variant. GATK _VariantsToTable_ produces a file that looks like: 
 
    CHROM    POS ID      QUAL    AC
     1        10  .       50      1
     1        20  rs10    99      10

Such a table contains only genomic (_i.e._ DNA-level) coordinates. Often the next question is what specific gene and peptide-level variants is each variant associated with? `cerebra` queries a reference genome (.fa) and annotation (.gtf) to match each DNA-level variant with its associated gene, probable transcript and probable peptide-level level variants. `cerebra` produces the following outfile: 

```
$ python
> import json
> f = open(/pth/to/cerebra/output.json)
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
    },
    "GOLGB1": {
        "A16_B000563": [],
        "A1_B001546": [
            "ENSP00000484083.1:p.?",
            "ENSP00000377275.3:p.(Ala1826Val)",
            "ENSP00000484083.1:p.(Gly1690Asp)",
            "ENSP00000484083.1:p.(Ala1746Val)",
            "ENSP00000377275.3:p.(Gly1770Asp)",
            "ENSP00000341848.5:p.(Thr911Ser)",
            "ENSP00000417767.1:p.(Thr782Ser)",
        ],
        "A1_B002531": [],
        "A1_B002570": [],
        "A2_B002558": [],
        "A3_B000561": [],
        "A3_B000568": [],
        "A3_B001544": [
            "ENSP00000377275.3:p.?",
            "ENSP00000341848.5:p.?",
            "ENSP00000484083.1:p.?"
        ],
        "A3_B002090": [],
        "A3_B002531": []
    }
}
```

Here _CCN1_ and _GOLGB1_ are gene names while _A16_B000563_, _A1_B001546_, _A1_B002531,_... are RNA-seq sample IDs. `cerebra` reports variants to every gene in the genome, for every sample in a given experiment. The _ENSP****_ numbers refer to [Ensembl](https://www.ensembl.org/index.html) translation IDs -- unique identifiers that correspond to exactly one polypeptide in the Ensembl database. The strings enclosed in parentheses refer to the amino-acid level variant reported in that particular sample. For example the string `Arg209Trp` indicates that position 209 of this particular polypeptide should contain an _Arg_, but the experimental sample instead contains a _Trp_. `cerebra` adheres to HGVS sequence variant [nomenclature](https://varnomen.hgvs.org/) in reporting amino-acid variants.

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
If you do not have access to "control" tissue, then proceed directly to `count-variants` or `find-peptide-variants`.

### `count-variants`
The `count-variants` module reports the raw variant counts for every gene across every sample.
The output is a CSV file that contains counts for each sample versus every gene in the genome. 

### `find-peptide-variants`
The `find-peptide-variants` module reports the peptide-level consequence of variants in the genome.
If working  with cancer samples, the user has the option to filter out all variants that are not found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database and are therefore unlikely to be pathogenic.
VCF records are converted to peptide-level variants, and then [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, 
in acordance to the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature. 
The output is a heirarchically ordered text file (CSV or JSON) that reports the the Ensemble protein ID and the gene associated with each variant, for each experimental sample. 

Variant callers are known to produce a great deal of false positives; the `--report-coverage` option is designed to give the user a greater degree of confidence in individual variant calls. If indicated this option will report raw counts for variant
and wildtype reads at each variant loci. We reason that variants with a high degree of read
support are less likely to be false positives. 

We should stress that `find-peptide-variants` does not *definitively* report peptide-level variants but rather the *likely*
set of peptide variants. 
Definitively reporting protein variants requires knowledge of alternate splicing -- this represents an open problem in scRNA-seq.
For example, if a read picks up a variant in exon 2 of _geneA_, we can report each of the potential spliceforms of _geneA_ that contain exon 2, but we **cannot** infer which of those particular spliceforms are actually present in our sample. 
Thus we report all possible spliceforms; determining the spliceform landscape of an individual cell from scRNA-seq is outside the scope of this project. 


Note that all three modules take advantage of multiprocessing.
Thus `cerebra` should scale better to high-memory machines with more cores, though it has been designed to run on everyday hardware. 

Installation
------------

From [PyPi](https://pypi.org/project/cerebra/)   
```pip install cerebra```

With [Docker](https://hub.docker.com/r/lincolnharris/cerebra)     
```docker pull lincolnharris/cerebra```                 
            
            
With git clone and the python standard library [venv](https://docs.python.org/3/library/venv.html) module
```
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
python3 -m venv cerebra-dev
source cerebra-dev/bin/activate
pip3 install -e . 
```

With git clone and [conda](https://docs.conda.io/en/latest/)
``` 
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install -e . 
```


Usage
-----

`cerebra` should now be installed as a commandline executable. 
`$ cerebra` should return help information

```
Usage: cerebra  <command>

  a tool for fast and accurate summarizing of variant calling format (VCF)
  files

Options:
  -h, --help  Show this message and exit.

Commands:
  germline-filter    filter out common SNPs/indels between control/germline samples and samples of interest
  count-variants    count total number of variants in each sample, and report on a per-gene basis
  find-peptide-variants  report peptide-level SNPs and indels in each sample, and associated coverage
```

Note that the `-h` command will display usage information for each of the three commands. 

An example workflow might look like this:   

**Step 1:**     
`cerebra germline-filter --processes 2 --control_path /path/to/control/vcfs --experimental_path /path/to/experimental/vcfs --metadata /path/to/metadata/file --outdir /path/to/filtered/vcfs`   

**Step 2:**     
`cerebra count-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --refgenome /path/to/genome/annotation --outfile /path/to/output/file /path/to/filtered/vcfs/*` 

**Step 3:**          
`cerebra find-peptide-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --annotation /path/to/genome/annotation --genomefa /path/to/genome/fasta --report_coverage 1 --output /path/to/output/file /path/to/filtered/vcfs/*`


Authors
--------
This work was produced by [Lincoln Harris](https://github.com/lincoln-harris), [Rohan Vanheusden](https://github.com/rvanheusden), [Olga Botvinnik](https://github.com/olgabot) and [Spyros Darmanis](https://spyrosdarmanis.wixsite.com/mylab) of the Chan Zuckerberg Biohub. For questions please contact ljharris018@gmail.com


Contributing
--------
We welcome any bug reports, feature requests or other contributions. Please submit a well documented report on our [issue tracker](https://github.com/czbiohub/cerebra/issues). For substantial changes please fork this repo and submit a pull request for review. 

See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/CONTRIBUTING.md) for additional details. 

You can find official releases [here](https://github.com/czbiohub/cerebra/releases). 
