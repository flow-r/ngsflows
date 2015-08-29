





## Downloading resource bundles.

Refernce libraries are available from Illumina and [Broad Insititute](https://www.broadinstitute.org/gatk/guide/article.php?id=1215).

We have a in built function to get genomes from Illumina (example below), but for Human genetics we recommend to the bundle from Broad Institute.

```
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


```
wget ... and download the reference data and files
## wget download the 
## test R
Rscript --version
## get runflow from github, name of the flow, samplesheet, parameterfile
## wget runflow
```


```r
install.packages('devtools')
require(devtools)
install_github(repo = 'ngsflows', username = 'sahilseth')
```
