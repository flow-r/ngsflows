

#' split_fq 
#' @param x file to split
#' @param split_size approx size of resulting file, in MB
#' @details
#' estimate the number of splits required
#' pre-compute the number of splits it would generate
#' make commands to make the split
#' approx, make 200mb files, for splitting...
split_fq <- function(x, split_size = 200){
  sz = file.size(x)/10^6
  
  ## if more than 1gb, only then splitting is worth it !
  if(split_size*2 > sz){
    message("No need to split")
    return()
  }
  chunks = round(sz/split_size)
  
  outbase = gsub("fastq.gz$|fastq$", "", x)
  splt_out = sprintf("%s%03d", outbase, 0:(tot_fls-1))
  final_out = paste0(splt_out, ".fastq")
}

## fqz=/rsrch2/iacs/ngs_runs/Y76I6Y76/Project_LTP9X-AMLLR/Sample_Futreal-LTP9X-AMLLR-240656/Futreal-LTP9X-AMLLR-240656_CTTGTA_L001_R1_001.fastq.gz
## perl ~/Dropbox/public/github_ngsflows/inst/scripts/fastq-splitter.pl --n-parts 278
##
##
##
##
##


if(FALSE){
  x = "/rsrch1/iacs/tmp/illumina_platinum/50x/NA12877/ERR194146.fastq.gz"
  #debug(split_fq)
  split_fq(x)
  
}


