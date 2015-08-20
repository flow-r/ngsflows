## --------- f l o w r        r e c i p i e s ------------- ##
## 
## 
## 

#' ln_bwa
#' @param fqs1 list of fastq files, may be a file with just the fastqs, one in each line.
#' @param fqs2 list of fastq files, may be a file with just the fastqs, one in each line. mate 2
#' @export
aln_bwa_merge <- function(fqs1, fqs2, 
  samplename,
  reflib = getOption("ngs_ref_bwa"),
  bwa_exe = getOption("ngs_bwa_exe"),
  bwa_aln_opts = getOption("ngs_bwa_aln_opts"),
  bwa_sampe_opts = getOption("ngs_bwa_sampe_opts"),
  samtools_exe = getOption("ngs_samtools_exe"),
  ## some platform
  seq_platform = "illumina", 
  rgpu = "lane1", center = "MDA"
  
  
  ){
    
    pipename = match.call()[[1]]
    message("Generating a ", pipename, " flowmat for sample: ", samplename)

    ## --- all subsequent steps would use this samplename
    ## --- code below looks a little cleaner
    options(ngs_samplename = samplename)
    
    ## --- Calling modules, each returns
    ##   - vector of outfiles
    ##   - a flowmat, which we need to rbind and are done !
    out_bwa = bwa(fqs1 = fqs1, fqs2 = fqs2)
    out_rg = picard_rg(out_bwa$outfiles)
    mergedbam = sprintf("%s.rg.sorted.bam", samplename)
    out_merge = picard_merge(out_rg$outfiles, mergedbam = mergedbam)

    ##--- creating flow mat
    flowmat = rbind(out_bwa$flowmat, out_rg$flowmat, out_merge$flowmat)

    return(list(outfile=mergedbam, flowmat = flowmat))
  }

## need 
## PE/SE
## add bowtie
## picard


if(FALSE){
  
  require(flowr)
  load_conf(search_conf("ngsflows_mda.conf"))
  
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




