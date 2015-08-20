## cpu_markdup=4, java_mem_markdup= "-Xmx8g",
## cpu_target=16,  java_mem_target= "-Xmx32g",
## cpu_realign=4,  java_mem_realign= "-Xmx4g", ## scatter 8 per node
## cpu_baserecalib=4,  java_mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
## cpu_printreads=4,  java_mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
#' Flow following Broad's best practices for variant calling, starting from sorted bam

#' @title Flow following Broad's best practices for variant calling, starting from sorted bam
#' @description This function provides a list of commands to run.
#' CPU's refer to number of cores reserved the job. We have a seperate param to control the number of threads the command would use. Example. On a node with 10GB and 10 cores. We may only reserve 10, but use only 2 threads for a 5GB memory command.
#' @param inbam
#' @param outbam
#' @param q_obj is provided output is a flow object, else a list of commands to run
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' @param cpu_markdup cpus to be reserved, default 4
#' @param java_mem_markdup memory provided to java
#' @param cpu_target
#' @param java_mem_target
#' @param cpu_realign
#' @param java_mem_realign
#' @param cpu_baserecalib
#' @param java_mem_baserecalib
#' @param cpu_printreads
#' @param java_mem_printreads
#' @param gatk_jar
#' @param picard_dir
#' @param reffa
#' @param gatk_target_opt
#' @param gatk_realign_opt
#' @param gatk_baserecalib_opt
#' @param printreads_opt
#' @export
bam_preprocess <- function(inbam, outbam, 
                           samplename = "samp",
                           split_by_chr = FALSE,
                           java_exe = get_opts("java_exe"),
                           java_tmp = get_opts("java_tmp"),
                           cpu_markdup = 4, mem_markdup= "-Xmx8g",
                           cpu_target = 16,  mem_target= "-Xmx32g",
                           cpu_realign = 8,  mem_realign= "-Xmx4g", ## scatter 8 per node
                           cpu_baserecalib = 4,  mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
                           cpu_printreads = 4,  mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
                           
                           gatk_jar = get_opts('gatk_jar'),
                           picard_dir = get_opts('picard_dir'),
                           
                           ref_fasta = get_opts('ref_fasta'),
                           
                           gatk_target_opts = get_opts('gatk_target_opts'),
                           gatk_realign_opts = get_opts('gatk_baserecalib_opts'),
                           gatk_baserecalib_opts = get_opts('gatk_baserecalib_opts'),
                           gatk_printreads_opts = get_opts('gatk_target_opts')){
  
  
  bam_prefix <- gsub(".bam", "", basename(bam))
  if(split_by_chr){
    chrs_info <- get_bam_chrs(bam)
    chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_")
  }else{
    chrs_prefix = bam_prefix
  }  
  
  ## ------------ dedup
  dedupbam <- paste0(bam_prefix, ".marked.bam")
  metricsfile <- paste0(bam_prefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/MarkDuplicates.jar INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                         java_exe, mem_markdup, java_tmp, picard_dir, bam, dedupbam, metricsfile)
  cmd_markdup
  
  ## ------------ realign
  intervalsfile <- paste0(bam_prefix, ".realign.intervals")
  realignedbams <- paste0(chrs_prefix ,".realigned.bam")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s",
                        java_exe, mem_target, java_tmp, gatk_jar, ref_fasta, dedupbam, 
                        intervalsfile, cpu_target, gatk_target_opts)
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar, ref_fasta, dedupbam, 
                         intervalsfile, realignedbams, gatk_realign_opts)
  
  ## ------------ base recalibration
  recalibbams <- gsub(".bam",".recalibed.bam",basename(inbam))
  recalibtabfile <- gsub(".bam",".recalib.tab",basename(inbam))
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar, ref_fasta, 
                             realignedbams, recalibtabfile,
                             gatk_baserecalib_opts)
  cmd_printreads <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s %s",
                            java_exe, mem_printreads, java_tmp, gatk_jar, ref_fasta, realignedbams, 
                            recalibtabfile, recalibbams,
                            gatk_printreads_opts)
  
  cmds <- list(markdup = cmd_markdup, 
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfile=recalibbams))
  
}

#' @importFrom Rsamtools scanBamHeader
#' @export
get_bam_chrs <- function(bam){
  out = scanBamHeader(bam)
  names(out[[1]]$targets)
  
}

.get_bam_chrs <- function(bam, samtools_exe = "samtools") {
  cmd <- sprintf("%s view -H %s | grep  '\\@SQ' | cut -f 2,3",
                 samtools_exe, bam)
  chrs_info <- system(cmd, intern = TRUE)
  chrs_info <- do.call(rbind, lapply(chrs_info, function(x){
    y = strsplit(x, "\t|:")[[1]]
    y[c(2,4)]
  }))
  return (chrs_info)
}


if(FALSE){
  bam="/scratch/iacs/gcc/levelii/140829_SN208_0523_BC5FKTACXX/TCGA-VD-AA8P-10A-01D-A40E_140829_SN208_0523_BC5FKTACXX_s_8_rg.sorted.bam"
  ##testing on
  ## /scratch/iacs/iacs_dep/sseth/flows/AdenoidCysticSanger/PD3176a/merge-preproc-0dd3b13a-aef5-4cdc-bee6-15db06f2dc71/tmp
}
