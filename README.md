
ngsflows
----------------------------------------------------

Provides a host of modules and pipelines supporting a host of genomic data ranging from 
whole genome/exome sequencing, RNA-Seq and chipseq etc.


```r
install.packages('devtools')
library(devtools)
install_github(repo = 'ngsflows', username = 'flow-r')
```


## Downloading resource bundles.

Reference libraries are available from Illumina and [Broad Insititute](https://www.broadinstitute.org/gatk/guide/article.php?id=1215).

We have a in built function to get genomes from Illumina (example below), but for Human genetics its better to use the bundle from Broad Institute.

One may run the following commands to get the bundle from Illumina


```
library(ngsflows)
fetch_genomes(species="Homo_sapiens")

################################################
Available Sources:
################################################

Ensembl
NCBI
UCSC

Example:
fetch_genomes species=Homo_sapiens src=NCBI build=build37.2

```

