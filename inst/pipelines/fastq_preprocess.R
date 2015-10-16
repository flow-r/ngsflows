
#' 
#' fastq_haplotyper
#' @param fqs1
#' @param fqs2
#' @param samplename
#' @param fqdir
#' 
#' @details 
#' This create a pipeline commands from fastq to variant calling
fastq_preprocess <- function(fqs1, fqs2, samplename = get_opts("samplename")){
  
  
  ## fetch the latest fastq_bam_bwa pipe
  #source("~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_bam_bwa.R")
  source(fetch_pipes('fastq_bam_bwa', silent = TRUE, last_only = TRUE)$pip)

  message("processing fastq_bam_bwa")
  f_merge = fastq_bam_bwa(fqs1, fqs2)
  message("processing bam_preprocess")
  f_preproc = preprocess(f_merge$outfile)

  flowmat = rbind(f_merge$flowmat, f_preproc$flowmat)
  
  outfiles = list(recalibed_bam = f_preproc$outfile)
  
  return(list(flowmat = flowmat, outfile = outfiles))
}


if(FALSE){
  
  ## load a configuration file with all the paths and options
  ## example of a cleaner approach
  set_opts(samplename = "samp1")
  load_opts(fetch_conf("ngsflows_mda.conf"), check = FALSE)
  get_opts("platform")
  
  ## creating a dummy flowmat
  out = fastq_preprocess(fqs1 = "my.fastq", fqs2 = "my.fastq")
  
  def = as.flowdef("inst/pipelines/fastq_preprocess.def")
  plot(def, pdffile = "inst/pipelines/fastq_preprocess.pdf")

  
}