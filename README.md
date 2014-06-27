---
title: "ngsutils manual"
author: "Sahil Seth"
date: "June 27, 2014"
output:
  pdf_document:
    toc: yes
---

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


### Setup
```
wget ... and download the reference data and files
## wget download the 
## test R
Rscript --version
## get runflow from github, name of the flow, samplesheet, parameterfile
## wget runflow
```


```r
#install.packages('devtools')
require(devtools)
```

```
## Loading required package: devtools
## 
## Attaching package: 'devtools'
## 
## The following objects are masked from 'package:utils':
## 
##     ?, help
## 
## The following object is masked from 'package:base':
## 
##     system.file
```

```r
install_github(repo = 'ngsutils', username = 'sahilseth')
```

```
## Installing github repo ngsutils/master from sahilseth
## Downloading master.zip from https://github.com/sahilseth/ngsutils/archive/master.zip
## Installing package from /var/folders/zk/n52wn3r928z88m0x66r5y8pny_mmbk/T//RtmpTmUi6B/master.zip
## arguments 'minimized' and 'invisible' are for Windows only
## Installing ngsutils
## '/Library/Frameworks/R.framework/Resources/bin/R' --vanilla CMD INSTALL  \
##   '/private/var/folders/zk/n52wn3r928z88m0x66r5y8pny_mmbk/T/RtmpTmUi6B/devtools379d625fa8e8/ngsutils-master'  \
##   --library='/Users/sseth/Rlibs' --install-tests
```
