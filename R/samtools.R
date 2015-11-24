#' @rdname samtools
#' @title set of wrapper functions around samtools
#' @param x input sam/bam file to be sorted and indexed
samtools_sort_index <- function(x, outbam, samtools_exe){

  cmd_sort = sprintf("%s sort %s -o %s",
    samtools_exe, x, gsub(outbam, ".bam", ""))

  cmd_index = sprintf("%s index %s",
    samtools_exe, outbam)

  cmds = list(sort = cmd_sort, index = cmd_index)
  flowmat = to_flowmat(cmds, samplename)
  return(list(flowmat = flowmat, outfiles = outbam) )

}
