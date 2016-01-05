

if(FALSE){

  require(ngsflows)
  fqdir = "/rsrch2/iacs/ngs_runs/150806_SN1120_0348_BC79KEACXX/fastqs/Y76I8n1Y76/Project_Pancreatic-MS132/Sample_GLizee-Pancreatic-MS132-MP013normalDNA"
  flowdef = as.flowdef("~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_haplotyper.def")
  plot_flow(flowdef)
  out = create_fq_sheet(x = fqdir)

  ## create a full flowmat

  fqs1 = head(subset(out, read == 1)$file, 2)
  fqs2 = head(subset(out, read == 2)$file, 2)
  ## get sample_name as a option for all subsequent steps
  opts_flow$set(samplename = "MS132")

  source('~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_haplotyper.R')
  load_opts('~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_haplotyper.conf', check = FALSE)
  out = fastq_haplotyper(fqs1, fqs2, samplename = s)

  fobj = to_flow(out$flowmat, flowdef, flowname = "fastq_haplotyper", flow_run_path = "flowr_test")
  fobj = submit_flow(fobj, execute = TRUE)

  tail(out$flowmat$cmd)
  fobj@jobs$printreads@cmds

  submit_flow(fobj, execute = TRUE)

  ## rerun for while debuggig
  wd = "/rsrch2/iacs/iacs_dep/sseth/flowr/runs/flowr_test/fastq_haplotyper-MS132-20150824-16-37-58-XScJT0OZ"
  fobj = rerun(x = wd, mat = out$flowmat, def = flowdef, start_from = "haplotyper", execute = TRUE)
  get_status(fobj)



}

#'
#' fastq_haplotyper
#' @param fqs1
#' @param fqs2
#' @param samplename
#' @param fqdir
#'
#' @details
#' This create a pipeline commands from fastq to variant calling
fastq_haplotyper <- function(fqs1, fqs2, samplename = opts_flow$get("samplename")){

  ## load a configuration file with all the paths and options
  ## example of a cleaner approach
  load_opts(fetch_conf("ngsflows_mda.conf"), check = FALSE)

  ## fetch the latest fastq_bam_bwa pipe
  #source("~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_bam_bwa.R")
  source(fetch_pipes('fastq_bam_bwa', silent = TRUE)$pip)

  message("processing fastq_bam_bwa")
  f_merge = fastq_bam_bwa(fqs1, fqs2)
  message("processing bam_preprocess")
  f_preproc = bam_preprocess(f_merge$outfile)
  message("processing haplotyper")
  f_hapl = haplotyper(f_preproc$outfile)

  flowmat = rbind(f_merge$flowmat, f_preproc$flowmat, f_hapl$flowmat)

  outfiles = list(recalibed_bam = f_preproc$outfile, variant_calls = f_hapl$outfile)

  return(list(flowmat = flowmat, outfile = outfiles))
}
