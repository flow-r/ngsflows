
# this pipelines imports the following:
# source(fetch_pipes('fastq_bam_bwa', last_only = TRUE)$pip)


#' fastq_haplotyper
#' @param fqs1
#' @param fqs2
#' @param samplename
#' @param fastq_path alternatively, one may supply a folder with fastq files
#'  for ONE sample. Alternatively, you may use \link{create_fq_sheet} to
#'  create a fastq sheet and loop over fastq_preprocess; for each sample.
#' @param fastq_format the fastq format to follow, see \link{create_fq_sheet}
#' and \link{split_names_fastq2} for more details
#' @param fqdir
#'
#' @details
#' This create a pipeline commands from fastq to variant calling
fastq_preprocess <- function(fqs1, fqs2,
                             fastq_path,
                             samplename = opts_flow$get("samplename"),
                             fastq_format = opts_flow$get("fastq_format")

){


  opts_flow$set(samplename = samplename)

  if(!missing(fastq_path)){
    # get a fastq sheet
    sheet = create_fq_sheet(x=fastq_path, format = fastq_format)
    sheet = subset(sheet, !is.na(num)) # remove rows, where num is NA. WHY??

    # extract list of fastq file from each
    if('read' %in% colnames(out)){

      fqs1 = subset(sheet, read == 1)$file
      fqs2 = subset(sheet, read == 2)$file
      lane = subset(sheet, read == 1)

      }else{
      fqs1 = out$file

      }

  }

  message("> processing fastq_bam_bwa ...")
  ## fetch the latest fastq_bam_bwa pipe
  #source("~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_bam_bwa.R")
  f_merge = fastq_bam_bwa(fqs1, fqs2)

  message("> processing bam_preprocess")
  f_preproc = preprocess(f_merge$outfiles)

  flowmat = rbind(f_merge$flowmat, f_preproc$flowmat)

  outfiles = list(recalibed_bam = f_preproc$outfile)

  return(list(flowmat = flowmat, outfiles = outfiles))
}


if(FALSE){

  ## load a configuration file with all the paths and options
  ## example of a cleaner approach
  opts_flow$load(fetch_conf("ngsflows_mda.conf"), check = FALSE)
  opts_flow$set(samplename = "samp1",
                rg_platform = "illumina",
                rg_center = "institute",
                rg_lane = 1,
                cpu_bwa_aln = 12)

  ## creating a dummy flowmat
  #debug(picard_rg)
  out = fastq_preprocess(fqs1 = "my.fastq", fqs2 = "my.fastq")


  def = as.flowdef("inst/pipelines/fastq_preprocess.def")
  plot(def, pdffile = "inst/pipelines/fastq_preprocess.pdf")


}
