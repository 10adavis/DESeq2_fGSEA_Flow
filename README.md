---
title: "README: Title of repo version 0.0"
author: "Andrew Davis"
date: "2021-12-17"
output:
  github_document:
    toc: TRUE
---

## Overview

This directory contains an analysis of....


The experiments were performed by ______ in the lab of ____.

The full analysis for this study can be found on_____

This analysis was version controlled as a GitHub repository, and can be accessed at the following url [](). Updates to this analysis will be tracked at the this location. For access to this repository, please reach out to ______.

## Lead contact(s):

* [Andrew J. Davis, PhD: 88adavis@gmail.com](mailto:88adavis@gmail.com)
* ....


## Metadata:

Metadata associated with this analysis can be found in the [Metadata](Metadata) folder. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Metadata/" from the [.gitignore](.gitignore) file.


## Input Data

The input data was acquired from ______ on  0/0/0000. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Input/" from the [.gitignore](.gitignore) file.

## Analysis pipeline:

1. 
2.
...


## Mechanics:

To run this analysis, first create/clean the results output folder by running the codes in [95_Make_Clean.Rmd](). Subsequently, run the code chunks in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd). This will run the Rmarkdown (.Rmd) files containing the actual code for this analysis in numerical order (i.e., 1.03_Gather-and-Tidy-Data.Rmd, followed by 2.03_Downstream_Analysis.Rmd). Note that one must modify the variable "files_in_r_to_run" in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) if one edits or add/deletes filenames of .Rmd scripts associated with this analysis. 

Running [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) will also render html files of each .Rmd file, which will be saved to the results folder, making useful reports of this analysis. Finally, this README.Rmd files will also be knitted to an html file, as well as a markdown (.md) file, in the working directory of this repository. This markdown file makes for easy viewing on GitHub, and acts as the "home page" for this repo.


## Output:

The resulting output files were saved to the [Results](Results) folder in this repository. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Results/" from the [.gitignore](.gitignore) file.
 

```
## [1] "1.03_Gather-and-Tidy-Data.html" "2.03_Downstream_Analysis.html"
```


## Summary: 

Describe any and all major insights generated from this analysis here....

Presentations and reports shared with other members of our team are stored in [Presentations_Reports](Presentations_Reports). Note that the contents of this folder will not be tracked by the remote repository by default (as these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Presentations_Reports/" from the [.gitignore](.gitignore) file.


## To do list:

1. 
...

## Template used:
This repository was generated from [10adavis/Rmarkdown_Template_2021](https://github.com/10adavis/Rmarkdown_Template).  

### Session information


```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19043)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] BiocManager_1.30.16 rmarkdown_2.11      here_1.0.1         
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.27     rprojroot_2.0.2   magrittr_2.0.1    evaluate_0.14    
##  [5] rlang_0.4.11      stringi_1.7.6     renv_0.14.0       jquerylib_0.1.4  
##  [9] tools_4.1.0       stringr_1.4.0     xfun_0.25         yaml_2.2.1       
## [13] compiler_4.1.0    htmltools_0.5.1.1 knitr_1.33
```

This document was processed on: 2021-12-17.




