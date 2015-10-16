

#' a wrapper around gatk haplotyper
#' @param x input bam
#' @export
haplotyper <- function(x, 
                       samplename = get_opts("samplename"),
                       outfile,
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
  
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)
  
  bam_prefix <- gsub(".bam", "", basename(x))
  if(split_by_chr){
    chrs_info <- get_fasta_chrs(ref_fasta)
    chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_")
    intervals_opts = paste0(" -L ", chrs_info)             ## interval files
    
  }else{
    chrs_prefix = paste0(bam_prefix, ".")
    intervals_opts = ""
    
  }  

  outvcf <- paste0(chrs_prefix, ".haplotyper.vcf")
  
  cmd_haplotyper <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T HaplotypeCaller -R %s -I %s -o %s -nct %s %s %s",
                        java_exe, java_mem, java_tmp, gatk_jar,
                        ref_fasta, x, outvcf, cpu_haplotyper, haplotyper_opts, intervals_opts)
  
  if(missing(outfile))
    outfile = paste0(bam_prefix, "_merged_haplotyper.vcf")
  
  cmd_merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T CatVariants -R %s -o %s -assumeSorted %s",
                      java_exe, java_mem, java_tmp, gatk_jar, ref_fasta, outfile, 
                      paste(" -V", outvcf, collapse = ""))

  cmds <- list(haplotyper = cmd_haplotyper, merge_haplotyper = cmd_merge)

  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles=list(all = outvcf, merged = outfile)))
}
