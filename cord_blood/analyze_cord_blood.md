HCA - Cord Blood Analysis
================

Packages and Shared Code
========================

``` r
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(DarkKinaseTools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BerginskiRMisc))
suppressPackageStartupMessages(library(here))
```

``` r
source(here('shared_functions.R'))
```

Data Download
=============

``` r
if (! file.exists(here('cord_blood','ica_cord_blood_h5.h5'))) {
  download.file('https://s3.amazonaws.com/preview-ica-expression-data/ica_cord_blood_h5.h5',
                here('cord_blood','ica_cord_blood_h5.h5'))
}
```

Data Loading
============

Loading this data takes some time and lots of memory (&gt; 16 Gigs).

``` r
cord_blood = read_h5_file_to_tidy('ica_cord_blood_h5.h5') %>%
  mutate(patient = gsub(".*(CB[[:digit:]]).*","\\1",barcode))
```

    ## Joining, by = "gene_index"

Exploritory Data Analysis
=========================

Let's get some quick summaries of the full data set, namely the total number of mapped reads per cell and the total number of unique genes per cell.

``` r
per_cell_summary  = cord_blood %>%
  group_by(barcode) %>%
  summarise(read_count = sum(counts),
            unique_genes = n())

ggplot(per_cell_summary,aes(x=read_count)) + 
  geom_histogram(binwidth=20) + 
  theme_berginski() +
  labs(x='Read Count Sum Per Cell',y='')
```

![](analyze_cord_blood_files/figure-markdown_github/summary_read_count-1.png) So their is a fair bit of skew in this data past the apparent double peak at the lower end of the distribution. Let's take a closer look at the lower end of the distribution:

``` r
ggplot(per_cell_summary,aes(x=read_count)) + 
  geom_histogram(binwidth=50) + 
  theme_berginski() +
  labs(x='Read Count Sum Per Cell',y='') +
  xlim(c(0,10000)) +
  geom_vline(xintercept = 800,color='red')
```

    ## Warning: Removed 2771 rows containing non-finite values (stat_bin).

![](analyze_cord_blood_files/figure-markdown_github/read_count_zoom-1.png)

There is a clear trough at 800 total mapped reads that I've marked out in red, that will be the lower threshold to toss out cells that didn't property adhere (maybe?) to the beads. The upper boundary doesn't seem as clear to me, but there is probably some barcode sets that represent multiple cells in a single droplet. How about a max read count of 7500?

What about the unique number of genes in each barcode:

``` r
ggplot(per_cell_summary,aes(x=unique_genes)) + 
  geom_histogram(binwidth=50) + 
  theme_berginski() +
  labs(x='Unique Genes Per Cell',y='') +
  geom_vline(xintercept = 275, color='red')
```

![](analyze_cord_blood_files/figure-markdown_github/gene_count-1.png)

I don't have a particularly compelling reason to filter on the number of unique genes found per cell, so I don't see any particular reason to filter on this value.

For good measure, let's also take a look at the binned scatterplot comparing the read counts and unique genes.

``` r
ggplot(per_cell_summary,aes(x=unique_genes,y=read_count)) + 
  geom_hex() + theme_berginski() +
  labs(x='Unique Gene Count',y="Read Count Sum")
```

![](analyze_cord_blood_files/figure-markdown_github/unnamed-chunk-2-1.png)

As expected, low unique gene counts somewhat line up with the read count sum.

Data Filtering and Results
==========================

Let's run the filtering to remove the problematic barcodes and then get out the kinase reads.

``` r
cord_blood_filt = cord_blood %>% filter_single_cell(read_total_range = c(800,7500))

cord_blood_kinases = cord_blood_filt %>% filter(ensembl_gene_id %in% all_kinases$ensembl_gene_id) %>%
  left_join(all_kinases)
```

    ## Joining, by = "ensembl_gene_id"

    ## Warning: Column `ensembl_gene_id` joining character vector and factor,
    ## coercing into character vector

``` r
cord_blood_dark_kinases = cord_blood_kinases %>% filter(class == "Dark")

per_kinase_stats = cord_blood_kinases %>% group_by(gene_names,class) %>%
  summarise(total_reads = sum(counts), 
            total_cells = length(unique(barcode)), 
            read_variance = var(counts))
```

Summary Results
---------------

Now let's run through some basic results:

-   505 unique kinases are represented in the data
-   0.71% of the mapped reads are to a kinase
-   100% of the bone marrow cells express at least one kinase
-   96% of the cells express a dark kinase

Summary Distributions
---------------------

``` r
ggplot(per_kinase_stats,aes(x=total_reads,y=stat(density),color=class, fill=class)) + 
  geom_histogram(position="dodge") + theme_berginski() +
  labs(color="Kinase\nType",fill="Kinase\nType",x="Total Reads",y="Density")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](analyze_cord_blood_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
ggplot(per_kinase_stats,aes(x=total_cells,y=stat(density),color=class, fill=class)) + 
  geom_histogram(position="dodge") + theme_berginski() +
  labs(color="Kinase\nType",fill="Kinase\nType",x="Total Cell Count",y="Density")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](analyze_cord_blood_files/figure-markdown_github/unnamed-chunk-4-1.png)

### Most Expressed Dark Kinases

| gene\_names | class |  total\_reads|  total\_cells|  read\_variance|
|:------------|:------|-------------:|-------------:|---------------:|
| STK17A      | Dark  |        266363|        137473|       1.6666167|
| STK17B      | Dark  |        162852|        112985|       0.5804515|
| PRPF4B      | Dark  |         63284|         53489|       0.2130445|
| RIOK3       | Dark  |         46033|         39190|       0.2320074|
| MKNK2       | Dark  |         45000|         38290|       0.2296777|
| CDK11A      | Dark  |         37004|         33874|       0.0976821|
| MARK3       | Dark  |         36241|         33400|       0.0899831|
| CDK12       | Dark  |         36004|         31602|       0.8299998|
| HIPK1       | Dark  |         35898|         32590|       0.1121307|
| CLK3        | Dark  |         33183|         30850|       0.0796969|

### Most Expressed Light Kinases

| gene\_names | class |  total\_reads|  total\_cells|  read\_variance|
|:------------|:------|-------------:|-------------:|---------------:|
| STK4        | Light |        182320|        123525|       0.6215469|
| LCK         | Light |        158178|         95457|       1.0458591|
| CLK1        | Light |        143878|        103528|       0.5732959|
| ATM         | Light |        131941|         95004|       0.4939435|
| JAK1        | Light |        115838|         86216|       0.4374471|
| CAMK4       | Light |        103991|         75278|       0.4714204|
| CSNK1A1     | Light |        102200|         81974|       0.2827944|
| AAK1        | Light |        101953|         78404|       0.3608784|
| BRD2        | Light |         76450|         63213|       0.2513626|
| ROCK1       | Light |         75514|         64002|       0.2011741|

Kinase Correlations
-------------------

I wrote a bit of code to do cross-correlation analysis on the expression levels of each of the kinases compared with one another, but I don't really know what to do with it. It's also time consuming, so I've commented it out.

``` r
# kinase_correlations = correlate_single_cell_read_counts(cord_blood_kinases)
# 
# saveRDS(kinase_correlations,here('cord_blood','cord_blood_kinases_correlations.rds'))
```
