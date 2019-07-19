cerebra
================================

[![image](https://img.shields.io/travis/%7B%7B%20cookiecutter.github_organization%20%7D%7D/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://travis-ci.org/%7B%7B%20cookiecutter.github_organization%20%7D%7D/%7B%7B%20cookiecutter.repo_name%20%7D%7D)


[![codecov](https://codecov.io/gh/%7B%7B%20cookiecutter.github_organization%20%7D%7D/%7B%7B%20cookiecutter.repo_name%20%7D%7D/branch/master/graph/badge.svg)](https://codecov.io/gh/%7B%7B%20cookiecutter.github_organization%20%7D%7D/%7B%7B%20cookiecutter.repo_name%20%7D%7D)

[![image](https://img.shields.io/pypi/v/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://pypi.python.org/pypi/%7B%7B%20cookiecutter.repo_name%20%7D%7D)


What is cerebra?
-------------------------------------

finds mutants in your scRNA-seq experiment

-   Free software: MIT license
-   Documentation: <https://>czbiohub.github.io/cerebra

Installation
------------

To install this code, clone this github repository and use pip to install

```
git clone <https://github.com/>czbiohub/cerebra.git 
cd cerebra 

# create a new conda environment
conda create -n cerebra python=3.6
conda activate cerebra

# The "." means "install *this*, the folder where I am now"
pip install -e . 
```

Usage
-----

cerebra should now be installed as a commandline executable. 
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

-   TODO

