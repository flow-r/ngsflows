## --------- f l o w r        r e c i p i e s ------------- ##
##
##
## --rg-id                        <string>    (read group ID)
## --rg-sample                    <string>    (sample ID)
## --rg-library                   <string>    (library ID)
## --rg-description               <string>    (descriptive string, no tabs allowed)
## --rg-platform-unit             <string>    (e.g Illumina lane ID)
## --rg-center                    <string>    (sequencing center name)
## --rg-date                      <string>    (ISO 8601 date of the sequencing run)
## --rg-platform                  <string>    (Sequencing platform descriptor)

## fetch fastq
## 
## 
##
## 
##


#' fastq_bam_bwa
#'
#' @param fqs1 list of fastq files, may be a file with just the fastqs, one in each line.
#' @param fqs2 list of fastq files, may be a file with just the fastqs, one in each line. mate 2
#'
#' @details
#' If fqs2 is missing, automatically use single end
#' @export
fastq_bam_bwa <- function(
  fqs1, fqs2,
  samplename){

    pipename = match.call()[[1]]
    message("Generating a ", pipename, " flowmat for sample: ", samplename)

    ## --- all subsequent steps would use this samplename
    set_opts(samplename = samplename)

    ## Calling modules, each returns
    ##   - a vector of outfiles
    ##   - a flowmat, which we need to rbind and are done !
    out_bwa = bwa(fqs1 = fqs1, fqs2 = fqs2)
    out_rg = picard_rg(out_bwa$outfiles)
    mergedbam = sprintf("%s.rg.sorted.bam", samplename) ## feel free to change this !
    out_merge = picard_merge(out_rg$outfiles, mergedbam = mergedbam)

    ##--- merging three flowmats
    flowmat = rbind(out_bwa$flowmat, out_rg$flowmat, out_merge$flowmat)

    return(list(outfile=mergedbam, flowmat = flowmat))
  }


if(FALSE){

  require(flowr)
  load_opts(search_conf("ngsflows_mda.conf"))

  ## This fails, extension seems weird
  flow_mat = aln_bwa_merge(fqs1 = rep("hello.fq.gz", 10),
    fqs2 = rep("hello.fq", 10),
    samplename = "smp")


  ## This fails, length is not the same, for paired end
  flow_mat = aln_bwa_merge(fqs1 = rep("hello.fq", 10),
    fqs2 = rep("hello.fq", 11),
    samplename = "smp")


  out_alnmerge = aln_bwa_merge(fqs1 = rep("hello.fq", 10),
    fqs2 = rep("hello.fq", 10),
    samplename = "smp")


}


main <- function(){
  out_alnmerge = aln_bwa_merge(
    fqs1 = rep("hello.fq", 10),
    fqs2 = rep("hello.fq", 10),
    samplename = "smp")

  flowmat =  out_alnmerge$flowmat
  def = to_flowdef(flowmat)
  plot_flow(to_flow(flowmat, def))

}


