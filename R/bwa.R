## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------



#' get extention of fastq files
get_fq_ext <- function(x){
  ext = gsub(".*(fastq.gz$|fastq$|fq$)", "\\1", x)
  return(ext)
}


#' @title  chk fq:
#' @description 
#' Files should have same extension. If paired, both should be of the same length
#' @param fqs1 a vector of file paths
#' @param fqs2 a vector of file paths
chk_fq <- function(fqs1, fqs2){
  paired_end = FALSE
  
  if(!missing(fqs2)){
    paired_end = TRUE
    ext = unique(get_fq_ext(c(fqs1, fqs2)))
  }else{
    ext = unique(get_fq_ext(fqs1))
  }
  
  if(length(ext) > 1)
    stop("fastq with different extenstions found, this seems troublesome !", ext)
  
  if(paired_end)
    if(length(fqs1) != length(fqs2))
      stop("Length of fqs1 and fqs2 dont match ! Exiting...\nfqs1\n:", 
        fqs1, "\n\nfqs2\n:", fqs2)
  
  cat_cmd = ifelse(ext == "gz", "zcat", "cat")
  
  return(list(ext = ext, paired_end = paired_end, cat_cmd = cat_cmd))
}




## --- define some default options, you may include them in ngsflows.conf
set_opts(
  bwa_aln_opts = "-l -k 2 -n 3",
  bwa_sampe_opts = "-o 1000",
  bwa_index = "bwa_index",
  bwa_samse_opts = "-o 1000"
)


#' Wrapper for BWA sequence alignment tool
#' 
#' @param fastq1 character vector with full paths to fastq files for mate 1.
#' @param fastq2 character vector
#' 
#' @return 
#' A list, with cmds as one of the variables which
#' contains system commands to be run.
#' 
#' @export
#' @importFrom flowr check_args
#' @details 
#' Default params picked up from ngsflow.conf
bwa <- function(
  fqs1, 
  fqs2,
  samplename = get_opts("samplename"),
  paired_end, ## auto detect it fastq2, is available
  
  bwa_exe = get_opts("bwa_exe"), 
  ref_bwa = get_opts("ref_bwa"),
  bwa_aln_opts = get_opts("bwa_aln_opts"),
  cpu_bwa_aln = get_opts("cpu_bwa_aln"),
  bwa_sampe_opts = get_opts("bwa_sampe_opts"),
  bwa_samse_opts = get_opts("bwa_samse_opts"),
  samtools_exe = get_opts("samtools_exe")
  
  ){
  
  
  ## --- some generic steps which may be done in case of 
  ## --- all FQ inputs !
  chkfq <- chk_fq(fqs1 = fqs1, fqs2 = fqs2)
  if(missing(paired_end))
    paired_end = chkfq$paired_end
  
  
  
  #source('~/Dropbox/public/github_flow/R/checkmat_assert.R')
  ## no arguments should be NULL
  check_args()

  ##  ------- Set up all the files which would be used.
  sai_files1 = file.path(gsub(chkfq$ext, "sai", basename(fqs1)))
  if(paired_end)
    sai_files2 = file.path( gsub(chkfq$ext, "sai", basename(fqs2)))
  
  ## These would be out files !. ALWAYS USE basename
  bam_files = file.path(gsub(chkfq$ext, "bam", basename(fqs1)))
  bam_prefix = gsub(".bam", "", bam_files)

  ## --- BWA aln and sampe
  #bwa_method <- match.arg(bwa_method)

  if(!paired_end){
    cmd_aln1 = sprintf("%s aln -t %s %s %s > %s", 
                       bwa_exe, cpu_bwa_aln, bwa_aln_opts, ref_bwa, fqs1, sai_files1)
    cmd_samse = sprintf("%s sampe %s %s %s %s| %s view -Shu - > %s",
      bwa_exe, bwa_samse_opts, ref_bwa, sai_files1, fqs1, samtools_exe, bam_prefix)
    cmds = list(aln1 = cmd_aln1, cmd_samse = cmd_samse)
  }
  
  if(paired_end){
    cmd_aln1 = sprintf("%s aln -t %s %s %s %s > %s", 
      bwa_exe, cpu_bwa_aln, bwa_aln_opts, ref_bwa, fqs1, sai_files1)
    cmd_aln2 = sprintf("%s aln -t %s %s %s %s > %s",
      bwa_exe, cpu_bwa_aln, bwa_aln_opts, ref_bwa, fqs2, sai_files2)
    cmd_sampe = sprintf("%s sampe %s %s %s %s %s %s | %s view -Shu - | %s sort - %s",
      bwa_exe, bwa_sampe_opts, ref_bwa, sai_files1, sai_files2, fqs1, fqs2, samtools_exe, samtools_exe, bam_prefix)
    ## --- make a named list of commands
    cmds = list(aln1 = cmd_aln1, aln2 = cmd_aln2, sampe = cmd_sampe)
  }
  
  
  ## --- convert to flow_mat compatible DF.
  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  ret = list(outfiles = bam_files, flowmat = flowmat)
  return(ret)
}


bwa_mem <- function(){
  ## ---- NOT implemented !
  
  if(bwa_method == "mem" & paired_end){
    cmds <- sprintf("%s mem %s %s %s %s",
      bwa_exe, bwa_mem_opt, bwa_index, fastq1, fastq2)
  }else if(bwa_method == "mem" & !paired_end){
    cmds <- sprintf("%s mem %s %s",
      bwa_exe, bwa_mem_opt, bwa_index, fastq1)
  }
  
  
  ## --- convert to flow_mat compatible DF.
  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  ret = list(outfiles = bam_files, flowmat = flowmat)
  return(ret)
  
}


# 
# #require(flow)
# setClass("bwa", contains = "job",
#   representation(fastq1 = "character", ## submit job
#     fastq2 = "character", ## type of queue
#     paired_end = "logical", 
#     bwa_exe = "character",
#     bwa_command = "character",
#     bwa_ref = "character",
#     bwa_opt = "character"
#   )) ## address of head node

if(FALSE){
  #bwa()@cmds
  
  ## PE
  #debug(bwa)
  out = bwa(fqs1 = rep("read1.fq", 20),
    fqs2 = rep("read2.fq", 20))
  
  
  ## SE
  out = bwa(fqs1 = rep("read1.fq", 20))
  
  #cmd_aln1 <- paste(bwapath, "/bwa aln ",bwa_aln_opts," ",reflib," ",fqs1,sep = "")
  
}