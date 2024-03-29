
cppd: Capture PBAT Probe Designer
=================================

<!-- badges: start -->
<!-- badges: end -->
The goal of cppd is to provide tools to design capture-PBAT probes. 

Installation
------------

You can install cppd by running:

``` r
remotes::install_github("tanaylab/cppd")
```

Example
-------

``` r
library(cppd)
library(tidyverse)

# Set misha database
gsetroot("/net/mraid14/export/data/db/tgdb/mm9/trackdb/")

# Define regions (e.g. promoters)
genes <- gintervals.load('intervs.global.tss') %>% sample_n(1e4) %>% pull(geneSymbol)
gene_promoters <- cppd.gene_promoters(gset$name) %>% select(chrom, start, end, strand, gene=geneSymbol)
fwrite(gene_promoters, "regions.csv")

# create the config file
cppd.dump_example_config("config.yaml")

# -- fill the yaml file  --
file.edit("config.yaml")

# generate probes
probes <- cppd.generate_probes("config.yaml", return_probes = TRUE)
```
