RNA-seq Analysis
================
Bob Policastro

``` r
library("tidyverse")
library("DESeq2")
library("ComplexHeatmap")
library("PCAtools")
library("EnhancedVolcano")
```

### Preparing Data

#### Loading Sample Info

``` r
runinfo <- "run_info.tsv" %>%
  read_tsv(col_names=c("info", "SRR_1", "SRR_2")) %>%
  separate(info, into=c("GSM", "name", "organism", "seqtype"), sep="[;:] ")

knitr::kable(runinfo)
```

| GSM        | name            | organism                | seqtype | SRR\_1     | SRR\_2     |
|:-----------|:----------------|:------------------------|:--------|:-----------|:-----------|
| GSM2817318 | orgR\_GO\_m\_r4 | Drosophila melanogaster | RNA-Seq | SRR6181119 | SRR7243937 |
| GSM2817317 | orgR\_GO\_f\_r4 | Drosophila melanogaster | RNA-Seq | SRR6181118 | SRR7243936 |
| GSM2647277 | orgR\_GO\_m\_r3 | Drosophila melanogaster | RNA-Seq | SRR7243507 | SRR5639586 |
| GSM2647276 | orgR\_GO\_m\_r2 | Drosophila melanogaster | RNA-Seq | SRR7243506 | SRR5639585 |
| GSM2647275 | orgR\_GO\_m\_r1 | Drosophila melanogaster | RNA-Seq | SRR7243505 | SRR5639584 |
| GSM2647274 | orgR\_GO\_f\_r3 | Drosophila melanogaster | RNA-Seq | SRR7243504 | SRR5639583 |
| GSM2647273 | orgR\_GO\_f\_r2 | Drosophila melanogaster | RNA-Seq | SRR7243503 | SRR5639582 |
| GSM2647272 | orgR\_GO\_f\_r1 | Drosophila melanogaster | RNA-Seq | SRR7243502 | SRR5639581 |

#### Loading Counts

``` r
counts <- "results/counts/counts.tsv" %>%
  read_tsv(skip=1) %>%
  select(!c(Chr, Start, End, Strand, Length)) %>%
  rename_with(~str_extract(.x, "(?<=aligned/)[[:alnum:]]+(?=_Aligned)"), !Geneid)

knitr::kable(counts[1:5, 1:5])
```

| Geneid         | SRR5639586 | SRR7243504 | SRR7243502 | SRR7243505 |
|:---------------|-----------:|-----------:|-----------:|-----------:|
| Myo81F         |         20 |          0 |          0 |          1 |
| CR41571        |          0 |          0 |          0 |          0 |
| CR12798        |          0 |          0 |          1 |          0 |
| lncRNA:CR46123 |          1 |         23 |         27 |          0 |
| lncRNA:CR46122 |         28 |          0 |          1 |          0 |

#### Prettier Sample Names

``` r
runinfo <- runinfo %>%
  pivot_longer(starts_with("SRR"), names_to="SRR", values_to="accession") %>%
  transmute(accession, name=str_c(accession, name, sep="_"))

counts <- counts %>%
  pivot_longer(!Geneid, names_to="accession", values_to="count") %>%
  right_join(runinfo, by="accession") %>%
  select(!accession) %>%
  pivot_wider(names_from=name, values_from=count)

knitr::kable(counts[1:5, 1:5])
```

| Geneid         | SRR5639586\_orgR\_GO\_m\_r3 | SRR7243504\_orgR\_GO\_f\_r3 | SRR7243502\_orgR\_GO\_f\_r1 | SRR7243505\_orgR\_GO\_m\_r1 |
|:---------------|----------------------------:|----------------------------:|----------------------------:|----------------------------:|
| Myo81F         |                          20 |                           0 |                           0 |                           1 |
| CR41571        |                           0 |                           0 |                           0 |                           0 |
| CR12798        |                           0 |                           0 |                           1 |                           0 |
| lncRNA:CR46123 |                           1 |                          23 |                          27 |                           0 |
| lncRNA:CR46122 |                          28 |                           0 |                           1 |                           0 |

### Differential Gene Expression

#### Create Count Matrix

``` r
count_matrix <- counts %>%
  column_to_rownames("Geneid") %>%
  as.matrix

knitr::kable(count_matrix[1:5, 1:5])
```

|                | SRR5639586\_orgR\_GO\_m\_r3 | SRR7243504\_orgR\_GO\_f\_r3 | SRR7243502\_orgR\_GO\_f\_r1 | SRR7243505\_orgR\_GO\_m\_r1 | SRR5639581\_orgR\_GO\_f\_r1 |
|:---------------|----------------------------:|----------------------------:|----------------------------:|----------------------------:|----------------------------:|
| Myo81F         |                          20 |                           0 |                           0 |                           1 |                           0 |
| CR41571        |                           0 |                           0 |                           0 |                           0 |                           0 |
| CR12798        |                           0 |                           0 |                           1 |                           0 |                           1 |
| lncRNA:CR46123 |                           1 |                          23 |                          27 |                           0 |                          27 |
| lncRNA:CR46122 |                          28 |                           0 |                           1 |                           0 |                           1 |

#### Create Design Table

``` r
coldata <- count_matrix %>%
  colnames %>%
  {tibble(sample=.)} %>%
  separate(sample, into=c("accession", "run"), sep="(?<=^SRR[[:digit:]]{7})_") %>%
  mutate(name=str_extract(run, "[[:alnum:]_]+(?=_r[[:digit:]]$)")) %>%
  mutate(name=factor(name, levels=c("orgR_GO_f", "orgR_GO_m"))) %>%
  group_by(run) %>%
  mutate(run=str_c(name, sprintf("run_%s", cur_group_id()), sep="_")) %>%
  column_to_rownames("accession")

knitr::kable(coldata)
```

|            | run                 | name        |
|:-----------|:--------------------|:------------|
| SRR5639586 | orgR\_GO\_m\_run\_7 | orgR\_GO\_m |
| SRR7243504 | orgR\_GO\_f\_run\_3 | orgR\_GO\_f |
| SRR7243502 | orgR\_GO\_f\_run\_1 | orgR\_GO\_f |
| SRR7243505 | orgR\_GO\_m\_run\_5 | orgR\_GO\_m |
| SRR5639581 | orgR\_GO\_f\_run\_1 | orgR\_GO\_f |
| SRR7243503 | orgR\_GO\_f\_run\_2 | orgR\_GO\_f |
| SRR6181119 | orgR\_GO\_m\_run\_8 | orgR\_GO\_m |
| SRR7243936 | orgR\_GO\_f\_run\_4 | orgR\_GO\_f |
| SRR7243506 | orgR\_GO\_m\_run\_6 | orgR\_GO\_m |
| SRR5639584 | orgR\_GO\_m\_run\_5 | orgR\_GO\_m |
| SRR5639582 | orgR\_GO\_f\_run\_2 | orgR\_GO\_f |
| SRR6181118 | orgR\_GO\_f\_run\_4 | orgR\_GO\_f |
| SRR7243937 | orgR\_GO\_m\_run\_8 | orgR\_GO\_m |
| SRR7243507 | orgR\_GO\_m\_run\_7 | orgR\_GO\_m |
| SRR5639585 | orgR\_GO\_m\_run\_6 | orgR\_GO\_m |
| SRR5639583 | orgR\_GO\_f\_run\_3 | orgR\_GO\_f |

#### Differential Expression

``` r
dds <- count_matrix %>%
  DESeqDataSetFromMatrix(coldata, design=~name) %>%
  collapseReplicates(coldata$run) %>%
  DESeq

degs <- dds %>%
  results %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  as_tibble %>%
  arrange(padj)

knitr::kable(head(degs))
```

| gene    |  baseMean | log2FoldChange |     lfcSE |      stat | pvalue | padj |
|:--------|----------:|---------------:|----------:|----------:|-------:|-----:|
| abs     |  1408.644 |       2.459230 | 0.0621057 |  39.59749 |      0 |    0 |
| Sfxn1-3 |  4431.279 |       3.193999 | 0.0438852 |  72.78075 |      0 |    0 |
| CG9855  |  4280.192 |       3.071651 | 0.0506254 |  60.67414 |      0 |    0 |
| srl     |  3420.496 |      -4.114378 | 0.0739988 | -55.60063 |      0 |    0 |
| CG31523 |  6252.971 |       2.072156 | 0.0398310 |  52.02372 |      0 |    0 |
| ctrip   | 31675.557 |       2.374883 | 0.0569352 |  41.71202 |      0 |    0 |

#### Normalized Counts

``` r
rlog_counts <- dds %>%
  rlog(blind=FALSE) %>%
  assay %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  as_tibble

knitr::kable(rlog_counts[1:5, 1:5])
```

| gene           | orgR\_GO\_f\_run\_1 | orgR\_GO\_f\_run\_2 | orgR\_GO\_f\_run\_3 | orgR\_GO\_f\_run\_4 |
|:---------------|--------------------:|--------------------:|--------------------:|--------------------:|
| Myo81F         |          -2.2746283 |           -2.174007 |          -0.0237437 |          -0.0292518 |
| CR41571        |           0.0000000 |            0.000000 |           0.0000000 |           0.0000000 |
| CR12798        |          -0.5278542 |           -2.315698 |          -2.2492158 |          -2.2518834 |
| lncRNA:CR46123 |           5.3051718 |            5.444921 |           5.5372580 |           5.5581812 |
| lncRNA:CR46122 |           0.6305740 |           -2.152461 |          -2.0230778 |          -2.0282693 |

### Plots

#### PCA

``` r
metadata <- rlog_counts %>%
  column_to_rownames("gene") %>%
  colnames %>%
  {data.frame(type=ifelse(str_detect(., "_m_"), "male", "female"))} %>%
  magrittr::set_rownames(colnames(column_to_rownames(rlog_counts, "gene")))
  
pca_results <- rlog_counts %>%
  column_to_rownames("gene") %>%
  pca(metadata=metadata)
```

PCA biplot.

``` r
biplot(pca_results, colby="type")
```

![](/N/slate/rpolicas/CoffeeAndCodingIUB/rnaseq_example/results/markdown/rnaseq_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

PCA loadings.

``` r
plotloadings(pca_results, rangeRetain=0.01)
```

![](/N/slate/rpolicas/CoffeeAndCodingIUB/rnaseq_example/results/markdown/rnaseq_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Volcano

``` r
EnhancedVolcano(degs, degs$gene, "log2FoldChange", "padj")
```

    ## Warning: One or more p-values is 0. Converting to 10^-1 * current lowest non-
    ## zero p-value...

![](/N/slate/rpolicas/CoffeeAndCodingIUB/rnaseq_example/results/markdown/rnaseq_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
