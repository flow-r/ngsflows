


#' @title
#' fastq_bam_bwa
#'
#' @param fqs1 list of fastq files, may be a file with just the fastqs, one in each line.
#' @param fqs2 list of fastq files, may be a file with just the fastqs, one in each line. mate 2
#'
#' @return returns a single merged bam file
#'
#' @details
#' If fqs2 is missing, automatically use single end
#'
#' @export
fastq_bam_bwa <- function(fqs1, fqs2,
                          outfile,
                          samplename = opts_flow$get("samplename")){


  ## --- all subsequent steps would use this samplename
  check_args(ignore = c("outfile", "fqs2"))
  opts_flow$set(samplename = samplename)
  pipename = "fastq_bam_bwa"
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  ## Calling modules, each returns
  ##   - a vector of outfiles
  ##   - a flowmat, which we need to rbind and are done !
  out_bwa = bwa.backtrack(fqs1 = fqs1, fqs2 = fqs2)
  out_rg = picard_rg(out_bwa$outfiles)

  if(missing(outfile))
    outfile = sprintf("%s.rg.sorted.bam", samplename) ## feel free to change this !

  out_merge = picard_merge(out_rg$outfiles, mergedbam = outfile)

  ##  merging three flowmats ---
  flowmat = rbind(out_bwa$flowmat, out_rg$flowmat, out_merge$flowmat)

  return(list(outfiles = outfile, flowmat = flowmat))
}

## ----------------------

if(FALSE){
	## example
  require(flowr)
  load_opts(fetch_conf("fastq_bam_bwa.conf"))

  ## This fails, extension seems weird
  flow_mat = fastq_bam_bwa(fqs1 = rep("hello.fq.gz", 10),
                           fqs2 = rep("hello.fq", 10),
                           samplename = "smp")


  ## This fails, length is not the same, for paired end
  flow_mat = fastq_bam_bwa(fqs1 = rep("hello.fq", 10),
                           fqs2 = rep("hello.fq", 11),
                           samplename = "smp")

  ## this works
  out = fastq_bam_bwa(fqs1 = rep("hello.fq", 10),
                      fqs2 = rep("hello.fq", 10),
                      samplename = "smp")

  debug(bwa.backtrack)
  out = fastq_bam_bwa(fqs1 = rep("hello.fq", 10),
  										samplename = "smp")



}
