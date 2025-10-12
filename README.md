---
title: "Enter analysis name here"
author: "Andrew J. Davis, Ph.D."
date: "21 February 2025"
output:
  rmarkdown::html_document:
    toc: TRUE
---



DESeq2_fGSEA_Flow v1.0

## Introduction
DESeq2_fGSEA_Flow has been produced to provide a standardized RNA-seq analysis pipeline.

This pipeline utilizes [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [fGSEA](https://bioconductor.org/packages/release/bioc/html/fgsea.html) R packages to generate the following outputs:

-   Normalized count data
-   Quality control including metrics and plots
-   Differentially expressed genes (with statistics) between treatments/conditions
-   Visualization of DE results as volcano plots
-   Gene set enrichment results using the fGSEA package

DESeq2_fGSEA_Flow was developed to automate these analyses in a logical and fully reproducible manner.

To support the reproducibility of these analyses, this pipeline employs literate programming using [Rmarkdown(Rmd)](https://rmarkdown.rstudio.com/) and full R package dependency management using the [renv](https://rstudio.github.io/renv/articles/renv.html) package.

![](Images/pipeline_workflow.png)

By default this pipeline will run on the example RNA-seq data in the [test_data](test_data) folder.

This data consists of RNA-seq data generated in the following study: [Interferon Receptor Signaling Pathways Regulating PD-L1 and PD-L2 Expression. A. Garcia-Diaz, et al. 2017. Cell Reports](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)30525-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124717305259%3Fshowall%3Dtrue#secsectitle0035) and deposited on GEO, here: [GSE96619](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96619).

It is expected that this pipeline will be run on user provided data, in an iterative manner, with the user modifying the Rmarkdown scripts in folder [R_Scripts](R_Scripts).To do so, the user needs to input 1) Raw count data and 2) Associated sample metadata. It is expected that the user has a thorough understanding of RNA-seq and DESeq2 [M. Love, W. Huber, S. Anders. Genome Biology. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).

## Reasons to use RNAseq_GSEA_Flow:

1.  Run differential gene expression analysis on RNA-seq data using DESeq2, which is one of the most commonly used methods in bioinformatics.
2.  Quickly understand biological changes in your data, by running fGSEA, a popular R packages that performs fast, pre-ranked gene set enrichment analysis.
3.  Easily and quickly upload your own RNA-seq count data and sample metadata (drag and drop to the Input folder).
4.  Reminds you to capture essential metadata, in a consistent format, in this readme file - simply fill in the blanks in the pre-drafted report below.
5.  Can be run on a typical PC, on most operating systems, including Windows, Mac, and Linux, for free.
6.  Developed and shared as a fully open source, version controlled repository using GitHub, the most popular code hosting platform for version control and collaboration.
7.  Analyses can be quickly and easily reproduced by any other user when used in conjunction with git and GitHub. 
8. Automatically export the results as a single file for quick upload to Ingenuity Pathway Analysis (IPA) as a single dataset.

## Requirements:

-   R version 4.3.1 (2023-06-16 ucrt)
-   [RStudio v1.4 or greater](https://www.rstudio.com/products/rstudio/)
-   [pandoc](https://pandoc.org/installing.html)

In addition, the user must install the following R packages:

-   rmarkdown_2.9
-   knitr_1.33

## Installation

This pipeline is deployed as a github repository template. To apply this pipeline to your particular experiment, one should first create a new repository using this template.

![](Images/template_repo.png)

After creating a new Rproj using this template, the first required step is to install the R dependencies required to run the analysis (written in the .Rmd files in R_scripts). To do so, execute the following command in your R environment.


```r
renv::restore()
```

This will automatically install the packages required to run this analysis. By doing this, one will be able to execute the analysis using the same R packages (and the same versions) each time. Note that this will only need to be run once. If you require additional packages, or need to update the associated packages, be sure to execute the following command to update the [renv.lock](renv.lock) file:


```r
renv::snapshot()
```

### Modifying input datasets and setting DESEq2 design:

By default this pipeline will run on the example RNA-seq data in the [Example_data](Example_data) folder.

In order to run this pipeline on your data please add your counts and metadata (both as .csv files) to the [Input](../Input) folder. For a more detailed description of the file formatting requirements, see the instructions in [1.03_Gather-and-Tidy-Data.Rmd](R_scripts/1.03_Gather-and-Tidy-Data.Rmd). One will have to modify the directory paths to the counts and metadata files in [1.03_Gather-and-Tidy-Data.Rmd](R_scripts/1.03_Gather-and-Tidy-Data.Rmd).

The user will also have to specify the design of the DESeq2 model (in [2.03-DESeq2.Rmd](R_scripts/2.03-DESeq2.Rmd)) which includes specifying the factors of interest. For example, the default model includes a simple multifactor design, where we measure the effect of the Treatment, controlling for the subject effect.


```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ subject + group)
```

The user will also have to specify the specific contrasts that they want to perform. By default the model runs 1 comparison: "anti_PDL_treated_tumor" vs "baseline_tumor"

To extract results for your particular dataset, modify/add these comparison\_ variables. Each comparison should be set to a variable (e.g., comparison_1), which equals a character vector, consisting of 3 character strings:

1.  The first string needs to correspond to the factor in design formula
2.  The name of the numerator level for the fold change
3.  The name of the denominator level for the fold change.


```r
comparison_1<-c("Treatment", "anti_PDL_treated_tumor", "baseline_tumor")
# The comparisons then need to be merged into a "list"
comparisons<-list(comparison_1)
```

### Capturing metadata:

Finally, the user is highly encouraged to fill in the details for their particular experiment/analysis, including at a minimum the information included below. The user is also encouraged to delete this introduction to the template/pipeline, and start the readme with the title of their analysis, followed by the Overview section below.

## Overview

This directory contains an analysis of....

The experiments were performed by \_\_\_\_\_\_ in the lab of \_\_\_\_.

The full analysis for this study can be found on\_\_\_\_\_

This analysis was version controlled as a GitHub repository, and can be accessed at the following url . Updates to this analysis will be tracked at the this location. For access to this repository, please reach out to \_\_\_\_\_\_.

## Lead contact(s):

-   [Andrew J. Davis, PhD: 88adavis\@gmail.com](mailto:88adavis@gmail.com)
-   ....

## Metadata:

Metadata associated with this analysis can be found in the [Metadata](Metadata) folder. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Metadata/" from the [.gitignore](.gitignore) file.

## Input Data

The input data was acquired from \_\_\_\_\_\_ on 0/0/0000. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Input/" from the [.gitignore](.gitignore) file.

## Analysis pipeline:

1.  In [1.03_Gather-and-Tidy-Data.Rmd](R_scripts/1.03_Gather-and-Tidy-Data.Rmd), one reads in the raw counts and metadata.
2.  In [2.03-DESeq2.Rmd](R_scripts/2.03-DESeq2.Rmd), one uses the DESeq2 R package to perform QC and exploratory analysis of the raw count data. One then uses DESeq2 to perform differential gene expression analysis.
3.  In [2.05-Volcano-Plots.Rmd](R_scripts/2.05-Volcano-Plots.Rmd), one generates volcano plots for each desired comparison.
4.  In [3.03-GSEA.Rmd](R_scripts/3.03-GSEA.Rmd), one performs Gene set enrichment analysis (GSEA) using the fgsea R package against gene pathways/lists from GO, Reactome, and Wikipathways for each of the comparisons.
5.  In [4.03-Export_to_IPA.Rmd](R_scripts/4.03-Export_to_IPA.Rmd), one tabulates the DE results for an overall summary of the number of DEGs. ONe also generates and exports a table that can be uploaded into IPA for a Core Analysis.
6.  In [5.03-Gene_Heatmaps.Rmd](R_scripts/5.03-Gene_Heatmaps.Rmd), one generates a heatmap of the users favorite genes. By default the script plots scaled expression of a few select genes.

## Mechanics:

To run this analysis, first create/clean the results output folder by running the codes in [95_Make_Clean.Rmd](). Subsequently, run the code chunks in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd). This will run the Rmarkdown (.Rmd) files containing the actual code for this analysis in numerical order (i.e., 1.03_Gather-and-Tidy-Data.Rmd, followed by 2.03_Downstream_Analysis.Rmd). Note that one must modify the variable "files_in_r\_to_run" in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) if one edits or add/deletes filenames of .Rmd scripts associated with this analysis.

Running [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) will also render html files of each .Rmd file, which will be saved to the results folder, making useful reports of this analysis. Finally, this README.Rmd files will also be knitted to an html file, as well as a markdown (.md) file, in the working directory of this repository. This markdown file makes for easy viewing on GitHub, and acts as the "home page" for this repo.

## Output:

The resulting output files were saved to the [Results](Results) folder in this repository. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Results/" from the [.gitignore](.gitignore) file.


```
#>  [1] "1.03_Gather-and-Tidy-Data.html"                                        
#>  [2] "2.03_DESeq2.html"                                                      
#>  [3] "2.07_Volcano-Plots.html"                                               
#>  [4] "3.03_GSEA.html"                                                        
#>  [5] "4.03_Export-to-IPA.html"                                               
#>  [6] "anti_PDL_treated_tumor_vs_baseline_tumor_DESeq2.csv"                   
#>  [7] "anti_PDL_treated_tumor_vs_baseline_tumor_pvalue_histogram.tiff"        
#>  [8] "anti_PDL_treated_tumor_vs_baseline_tumor_volcano_padj-0.1_log2FC-1.pdf"
#>  [9] "DEG_Summary_Table.csv"                                                 
#> [10] "DESeq2_dispersions_plot.pdf"                                           
#> [11] "Filtered_Normalized_Counts.csv"                                        
#> [12] "Filtered_Raw_Counts.csv"                                               
#> [13] "GSEA_anti_PDL_treated_tumor_vs_baseline_tumor.csv"                     
#> [14] "GSEA_Top-20-Gene-Sets_anti_PDL_treated_tumor_vs_baseline_tumor.pdf"    
#> [15] "Heatmap_Select_Genes.pdf"                                              
#> [16] "PCA_PC1-vs-PC2_rlog.tiff"                                              
#> [17] "PCA_PC1-vs-PC2_rlog_ellipses.tiff"                                     
#> [18] "PCA_PC1-vs-PC2_VST.tiff"                                               
#> [19] "PCA_PC1-vs-PC2_VST_ellipses.tiff"                                      
#> [20] "Sample_Correlation_rlog.pdf"                                           
#> [21] "Sample_Correlation_VST.pdf"
```

## IPA:

Results from this experiment ([IPA/Merged_Treatment_Comparisons_DESeq2_Results.csv](IPA/Merged_Treatment_Comparisons_DESeq2_Results.csv)) were uploaded into IPA, under project "WWWWWW". The Dataset is titled: "XXXXXX" and the Core Analysis is titled: "ZZZZZ".

The resulting output files were saved to the [IPA](IPA) folder in this repository. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/IPA/" from the [.gitignore](.gitignore) file.


```
#> [1] "Merged_Treatment_Comparisons_DESeq2_Results.csv"
```

Describe any and all major insights generated from this analysis here....

Presentations and reports shared with other members of our team are stored in [Presentations_Reports](Presentations_Reports). Note that the contents of this folder will not be tracked by the remote repository by default (as these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Presentations_Reports/" from the [.gitignore](.gitignore) file.

## To do list:

1.  ...

## Template used:

This repository was generated from [10adavis/DESeq2_fGSEA_Flow v1.0](https://github.com/10adavis/DESeq2_fGSEA_Flow).

### Session information


```r
sessionInfo()
#> R version 4.3.1 (2023-06-16 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19045)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices datasets  utils    
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] genefilter_1.82.1           gplots_3.1.3               
#>  [3] scales_1.2.1                msigdbr_7.5.1              
#>  [5] biomartr_1.0.5              data.table_1.14.8          
#>  [7] GSEABase_1.62.0             graph_1.78.0               
#>  [9] annotate_1.78.0             XML_3.99-0.14              
#> [11] reactome.db_1.84.0          GO.db_3.17.0               
#> [13] fgsea_1.26.0                dplyr_1.1.3                
#> [15] magrittr_2.0.3              EnhancedVolcano_1.18.0     
#> [17] ggrepel_0.9.4               rlist_0.4.6.2              
#> [19] pheatmap_1.0.12             org.Hs.eg.db_3.17.0        
#> [21] AnnotationDbi_1.62.2        readxl_1.4.3               
#> [23] ggplot2_3.4.4               ashr_2.2-63                
#> [25] DESeq2_1.40.2               SummarizedExperiment_1.30.2
#> [27] Biobase_2.60.0              MatrixGenerics_1.12.3      
#> [29] matrixStats_1.0.0           GenomicRanges_1.52.1       
#> [31] GenomeInfoDb_1.36.4         IRanges_2.34.1             
#> [33] S4Vectors_0.38.2            BiocGenerics_0.46.0        
#> [35] future_1.33.0               rmarkdown_2.29             
#> [37] here_1.0.1                 
#> 
#> loaded via a namespace (and not attached):
#>  [1] RColorBrewer_1.1-3      jsonlite_1.8.7          farver_2.1.1           
#>  [4] zlibbioc_1.46.0         vctrs_0.6.4             memoise_2.0.1          
#>  [7] RCurl_1.98-1.12         SQUAREM_2021.1          mixsqp_0.3-48          
#> [10] htmltools_0.5.6.1       S4Arrays_1.0.6          progress_1.2.2         
#> [13] curl_5.1.0              truncnorm_1.0-9         cellranger_1.1.0       
#> [16] sass_0.4.7              parallelly_1.36.0       KernSmooth_2.23-21     
#> [19] bslib_0.5.1             cachem_1.0.8            lifecycle_1.0.3        
#> [22] pkgconfig_2.0.3         Matrix_1.5-4.1          R6_2.5.1               
#> [25] fastmap_1.1.1           GenomeInfoDbData_1.2.10 digest_0.6.33          
#> [28] colorspace_2.1-0        rprojroot_2.0.3         irlba_2.3.5.1          
#> [31] RSQLite_2.3.1           invgamma_1.1            filelock_1.0.2         
#> [34] labeling_0.4.3          fansi_1.0.5             httr_1.4.7             
#> [37] abind_1.4-5             compiler_4.3.1          bit64_4.0.5            
#> [40] withr_2.5.1             BiocParallel_1.34.2     DBI_1.1.3              
#> [43] biomaRt_2.56.1          MASS_7.3-60             rappdirs_0.3.3         
#> [46] DelayedArray_0.26.7     gtools_3.9.4            caTools_1.18.2         
#> [49] tools_4.3.1             glue_1.6.2              grid_4.3.1             
#> [52] generics_0.1.3          snow_0.4-4              gtable_0.3.4           
#> [55] hms_1.1.3               xml2_1.3.5              utf8_1.2.3             
#> [58] XVector_0.40.0          pillar_1.9.0            stringr_1.5.0          
#> [61] babelgene_22.9          splines_4.3.1           BiocFileCache_2.8.0    
#> [64] lattice_0.21-8          survival_3.5-5          renv_1.0.3             
#> [67] bit_4.0.5               tidyselect_1.2.0        locfit_1.5-9.8         
#> [70] Biostrings_2.68.1       knitr_1.44              xfun_0.40              
#> [73] stringi_1.7.12          yaml_2.3.7              evaluate_0.22          
#> [76] codetools_0.2-19        tibble_3.2.1            cli_3.6.1              
#> [79] xtable_1.8-4            munsell_0.5.0           jquerylib_0.1.4        
#> [82] Rcpp_1.0.11             globals_0.16.2          dbplyr_2.3.4           
#> [85] png_0.1-8               blob_1.2.4              prettyunits_1.2.0      
#> [88] bitops_1.0-7            listenv_0.9.0           crayon_1.5.2           
#> [91] rlang_1.1.1             cowplot_1.1.1           fastmatch_1.1-4        
#> [94] KEGGREST_1.40.1
```

This document was processed on: 2025-02-21.
