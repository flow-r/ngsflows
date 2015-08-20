

#' @param fq1 a single fq file. Can be fastq.gz
#' @param fq2 a single fq file. Can be fastq.gz
#' @param chunks number of chunks into which both these files should be split
split_aln_merge <- function(fq1, fq2, chunks = 100, samplename){

  ## --- needs this other pipeline
  source(fetch_pipes("aln_bwa_merge")$pipe)

  ## to handle both reads, might need to handle the jobname
  out_split1 = split_fq(fq1, samplename = samplename,
  											chunks = chunks, jobname = "split_fq1")
  out_split2 = split_fq(fq2, samplename = samplename,
  											chunks = chunks, jobname = "split_fq2")

  ## supply both files, since this is paired end.
  ## use only one for single end
  out_aln_merge = aln_bwa_merge(fqs1 = out_split1$outfiles,
    fqs2 = out_split1$outfiles,
    samplename = samplename)

  finaloutfile = out_aln_merge$outfiles

  ## --- stitch all the commands together!
  flowmat = rbind(
    out_split1$flowmat,
    out_split2$flowmat,
    out_aln_merge$flowmat)

  return(list(flowmat = flowmat, outfiles = finaloutfile))
}


## --- Following runs a small example
main <- function(){
  require(flowr);require(ngsflows)
  pip_name = "split_aln_merge"

  fq1 = "/rsrch1/iacs/tmp/illumina_platinum/50x/NA12877/ERR194146_1.fastq.gz"
  fq2 = "/rsrch1/iacs/tmp/illumina_platinum/50x/NA12877/ERR194146_2.fastq.gz"
  samplename = "ERR194146"

  out = split_aln_merge(fq1 = fq1, fq2 = fq2)
  #flowdef = to_flowdef(out$flowmat)
  plot_flow(to_flow(out$flowmat, flowdef))
  def = as.flowdef(file.path(pip_path, "split_aln_merge.def"))

  fobj = to_flow(out$flowmat, def, flowname = pip_name)

  submit_flow(fobj)


}



