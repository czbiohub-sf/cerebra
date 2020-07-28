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
