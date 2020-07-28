After following the installation procedure outlined in INSTALL.md:  

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
