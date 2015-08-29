

#' a wrapper around gatk haplotyper
#' @param x input bam
#' @param x input bam
#' @param x input bam
#' @param x input bam
#' @param x input bam
#' @export
haplotyper <- function(x, 
                       samplename = get_opts("samplename"),
                       split_by_chr = TRUE,
                       
                       ref_fasta = get_opts("ref_fasta"),
                       java_exe = get_opts("java_exe"),
                       java_mem = get_opts("java_mem"),
                       java_tmp = get_opts("java_tmp"),
                       gatk_jar = get_opts("gatk_jar"),
                       haplotyper_opts = get_opts("haplotyper_opts"),
                       cpu_haplotyper = get_opts("cpu_haplotyper")) {
  
  ## no args should be null
  check_args()  
  
  bam_prefix <- gsub(".bam", "", basename(x))
  if(split_by_chr){
    chrs_info <- get_bam_chrs(x)
    chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_")
  }else{
    chrs_prefix = bam_prefix
  }  

  outvcf <- paste0(chrs_prefix, ".haplotyper.vcf")
  
  cmd_haplotyper <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T HaplotypeCaller -R %s -I %s -o %s -nct %s",
                        java_exe, java_mem, java_tmp, gatk_jar,
                        ref_fasta, x, outvcf, cpu_haplotyper)

  cmds <- list(haplotyper = cmd_haplotyper)

  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfile=outvcf))
}
