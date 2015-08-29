## cpu_markdup=4, java_mem_markdup= "-Xmx8g",
## cpu_target=16,  java_mem_target= "-Xmx32g",
## cpu_realign=4,  java_mem_realign= "-Xmx4g", ## scatter 8 per node
## cpu_baserecalib=4,  java_mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
## cpu_printreads=4,  java_mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
#' Flow following Broad's best practices for variant calling, starting from sorted bam

#' @title Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file
#' @description This function provides a wrapper around the best practices described on \href{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}{GATK's website}.
#' If the link is broken google 'GATK best practices'
#' 
#' This aims to perform the following steps ( for DNA ):
#' 
#' \itemize{
#' \item mark duplicates
#' \item realign indels
#' \item recalibrate bases
#' \item current version: \emph{3.4-46}
#' }
#' 
#' For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).
#' 
#' \strong{NOTE}:
#' 
#' Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{CPU threads while others use data threads},
#' flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
#' 
#' @param inbam
#' @param outbam
#' @param q_obj is provided output is a flow object, else a list of commands to run
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' @param cpu_markdup not used.
#' @param java_mem_markdup memory provided to java
#' @param cpu_target number of threads used for GATK target creation step
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
#' 
#' @examples \donotrun{
#' ## load options, including paths to tools and other parameters
#' load_opts(fetch_conf("ngsflows.conf"), check = FALSE)
#' out = bam_preprocess("my_wex.bam")
#' 
#' }
bam_preprocess <- function(x, 
                           outfile, 
                           samplename = get_opts("samplename"),
                           split_by_chr = FALSE,
                           java_exe = get_opts("java_exe"),
                           java_tmp = get_opts("java_tmp"),
                           cpu_markdup = 1,
                           mem_markdup= "-Xmx8g",
                           cpu_target = get_opts("cpu_target"),  ## not used
                           mem_target= "-Xmx32g",
                           cpu_realign = get_opts("cpu_realign"),
                           mem_realign= "-Xmx4g", ## scatter 8 per node
                           cpu_baserecalib = get_opts("cpu_baserecalib"),  
                           mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
                           cpu_printreads = get_opts("cpu_printreads"),
                           mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
                           
                           gatk_jar = get_opts('gatk_jar'),
                           picard_dir = get_opts('picard_dir'),
                           
                           ref_fasta = get_opts('ref_fasta'),
                           
                           gatk_target_opts = get_opts('gatk_target_opts'),
                           gatk_realign_opts = get_opts('gatk_realign_opts'),
                           gatk_baserecalib_opts = get_opts('gatk_baserecalib_opts'),
                           gatk_printreads_opts = get_opts('gatk_printreads_opts')){
  
  
  ## determine output file name
  if(missing(outfile))
    bam_prefix <- gsub(".bam", "", basename(x))
  else
    bam_prefix <- gsub(".bam", "", basename(outfile))
  
  ## if file is available determine whether to split for faster processing
  if(split_by_chr){
    #chrs_info <- get_bam_chrs(x)
    chrs_info <- get_fasta_chrs(ref_fasta)
    chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_") ## bam names
    intervals_opts = paste0(" -L ", chrs_info)             ## interval files
  }else{
    chrs_prefix = bam_prefix
    intervals_opts = ""
  }  
  
  check_args(ignore = "outfile")
  
  ## ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bam_prefix, ".marked.bam")
  metricsfile <- paste0(bam_prefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                         java_exe, mem_markdup, java_tmp, picard_dir, x, dedupbam, metricsfile)
  cmd_markdup
  
  ## ------------ realign; SINGLE FILE
  intervalsfiles <- paste0(chrs_prefix, ".realign.intervals")
  realignedbams <- paste0(chrs_prefix ,".realigned.bam")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s",
                        java_exe, mem_target, java_tmp, gatk_jar, ref_fasta, dedupbam, 
                        intervalsfiles, cpu_target, gatk_target_opts)
  
  
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar,          ref_fasta, dedupbam, 
                         intervalsfiles, realignedbams, gatk_realign_opts, intervals_opts)
  
  ## ------------ base recalibration
  recalibbams <- paste0(chrs_prefix, ".recalibed.bam")
  recalibtabfile <- paste0(chrs_prefix, ".recalib.tab")
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s -nct %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar, ref_fasta, 
                             realignedbams, recalibtabfile, cpu_baserecalib,
                             gatk_baserecalib_opts)
  cmd_printreads <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s -nct %s %s",
                            java_exe, mem_printreads, java_tmp, gatk_jar, ref_fasta, realignedbams, 
                            recalibtabfile, recalibbams, cpu_printreads,
                            gatk_printreads_opts)
  
  cmds <- list(markdup = cmd_markdup, 
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfile=recalibbams))
  
}

#' @export
get_bam_chrs <- function(x){
  if(file.exists(x)){
    out = Rsamtools:::scanBamHeader(x)
    chrs = names(out[[1]]$targets)
  }else{
    message("bam does not exists, returning hg19 chrs")
    chrs = c(1:22,"X","Y","MT")
  }
  return(chrs)
}


#' Read the associated dictionary file and return a list of chromosome names
#'
#' @param x a reference genome fasta file
#'
#' @export
#' 
#' @importFrom params read_sheet
#'
get_fasta_chrs <- function(x){
  dict = gsub("fasta$", "dict", x)
  if(!file.exists(dict))
    stop(c("We need a .dict for the reference fasta file to proceed.", 
         "Follow this link to learn more: http://lmgtfy.com/?q=create+dict+fasta"))
    seqs = read_sheet(dict, ext = "tsv", skip = 1)
    gsub("SN:", "", seqs[, 2], fixed = TRUE)
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
