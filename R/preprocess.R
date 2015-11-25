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
#' @param q_obj is provided output is a flow object, else a list of commands to run
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' @param cpu_markdup not used.
#' @param java_mem_markdup memory provided to java
#' @param cpu_target number of threads used for GATK target creation step
#'
#' @export
#'
#' @examples \dontrun{
#' ## load options, including paths to tools and other parameters
#' load_opts(fetch_conf("ngsflows.conf"), check = FALSE)
#' out = bam_preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)
#'
#' }
#'
preprocess <- function(x,
                       outfile,
                       samplename = opts_flow$get("samplename"),
                       split_by_chr = opts_flow$get("split_by_chr"),

                       java_exe = opts_flow$get("java_exe"),
                       java_tmp = opts_flow$get("java_tmp"),

                       gatk_jar = opts_flow$get('gatk_jar'),
                       picard_dir = opts_flow$get('picard_dir'),
                       samtools_exe = opts_flow$get('samtools_exe'),

                       cpu_markdup = 1,
                       mem_markdup= "-Xmx8g",
                       cpu_target = opts_flow$get("cpu_target"),  ## not used
                       mem_target= "-Xmx32g",
                       cpu_realign = opts_flow$get("cpu_realign"),
                       mem_realign= "-Xmx4g", ## scatter 8 per node
                       cpu_baserecalib = opts_flow$get("cpu_baserecalib"),
                       mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
                       cpu_printreads = opts_flow$get("cpu_printreads"),
                       mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8

                       ref_fasta = opts_flow$get('ref_fasta'),

                       gatk_target_opts = opts_flow$get('gatk_target_opts'),
                       gatk_realign_opts = opts_flow$get('gatk_realign_opts'),
                       gatk_baserecalib_opts = opts_flow$get('gatk_baserecalib_opts'),
                       gatk_printreads_opts = opts_flow$get('gatk_printreads_opts')){

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

  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  # ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bam_prefix, ".marked.bam")
  metricsfile <- paste0(bam_prefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                         java_exe, mem_markdup, java_tmp, picard_dir, x, dedupbam, metricsfile)
  cmd_markdup

  ## ------------ realign; SINGLE FILE
  intervalsfiles <- paste0(bam_prefix, ".realign.intervals")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s %s",
                        java_exe, mem_target, java_tmp, gatk_jar, ref_fasta, dedupbam,
                        intervalsfiles, cpu_target, gatk_target_opts)

  realignedbams <- paste0(chrs_prefix ,".realigned.bam")
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar, ref_fasta, dedupbam,
                         intervalsfiles, realignedbams, gatk_realign_opts, intervals_opts)

  ## ------------ base recalibration
  ## explicity define intervals here as well, though not required.
  recalibbams <- paste0(chrs_prefix, ".recalibed.bam")
  recalibtabfile <- paste0(chrs_prefix, ".recalib.tab")
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s -nct %s %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar, ref_fasta,
                             realignedbams, recalibtabfile, cpu_baserecalib,
                             gatk_baserecalib_opts, intervals_opts)
  cmd_printreads1 <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s -nct %s %s %s",
                            java_exe, mem_printreads, java_tmp, gatk_jar, ref_fasta, realignedbams,
                            recalibtabfile, recalibbams, cpu_printreads,
                            gatk_printreads_opts, intervals_opts)
  cmd_printreads2 <- sprintf("%s index %s", samtools_exe, recalibbams)
  cmd_printreads = sprintf("%s;%s", cmd_printreads1, cmd_printreads2)

  cmds <- list(markdup = cmd_markdup,
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)

  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles = recalibbams))

}

#' @export
get_bam_chrs <- function(x){
  if(file.exists(x)){
    out = Rsamtools::scanBamHeader(x)
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
#' @param length provide length of each chromosome as well [FALSE]
#'
#' @export
#'
#' @importFrom params read_sheet
#' @importFrom tools file_path_sans_ext
#'
get_fasta_chrs <- function(x = opts_flow$get("ref_fasta"),
                           length = FALSE){
  check_args()
  #dict = gsub("fasta$", "dict", x)
  dict = paste0(file_path_sans_ext(x), ".dict")
  if(!file.exists(dict)){
    message(c("We need a dictionary (extension: .dict) for the reference fasta file to proceed. ",
              "Follow this link to learn more: http://lmgtfy.com/?q=create+dict+fasta"))
    warning("By default using hg19 chrs.")
    chrs = c(1:22,"X","Y","MT")
  }else{
    seqs = read_sheet(dict, ext = "tsv", skip = 1, head = FALSE)
    chrs = gsub("SN:", "", seqs[, 2], fixed = TRUE)
    lens = gsub("LN:", "", seqs[, 3], fixed = TRUE)
    if(length)
      return(list(chrs = chrs, lens = lens))
  }
  return(chrs)
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
