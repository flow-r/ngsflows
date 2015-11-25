
#' A wrapper around somatic mutation caller MuTect
#'
#' This generates a set of commandlines, per chromosome
#'
#' @param tumor_bam path to a tumor bam file
#' @param normal_bam path to a normal bam file
#' @param samplename name of the sample, to be added to the flowmat
#' @param outfile output file name [optional]
#' @param is_merged specify if the input bam is already split by chromosomes (FALSE) or is a merged file with all chromosome (TRUE). [FALSE]
#' @param split_by_chr fork running mutect by chromosome to make it faster.
#' @param java_exe path to java's executable [java]
#' @param java_tmp path to a temp folder to be used by java
#' @param mutect_jar path to mutect's jar file
#' @param cpu_mutect integer specifying number of thread to be used per mutect fork
#' @param mem_mutect amount of memory to be used by mutect [-Xmx12g]
#' @param ref_fasta path to the reference genome in fasta format
#' @param mutect_opts additional arguments passed onto mutect
#'
#' @export
#'
#' @examples \dontrun{
#'
#' x = "tumor.bam"
#' y = "normal.bam"
#'
#' out = mutect(x, y, is_merged = TRUE)
#'
#' }
mutect <- function(tumor_bam, normal_bam,
                   samplename = opts_flow$get("samplename"),
                   outfile,

                   is_merged = TRUE,
                   split_by_chr = TRUE,

                   java_exe = opts_flow$get("java_exe"),
                   java_tmp = opts_flow$get("java_tmp"),

                   mutect_jar = opts_flow$get('mutect_jar'),

                   cpu_mutect = opts_flow$get('cpu_mutect'), ## curently not suported
                   mem_mutect = opts_flow$get("java_mem"),

                   ref_fasta = opts_flow$get('ref_fasta'),

                   mutect_opts = opts_flow$get('mutect_opts')

){

  ## determine output file name
  if(missing(outfile))
    bam_prefix <- gsub(".bam", "", basename(normal_bam))
  else
    bam_prefix <- outfile

    if(missing(is_merged))
      is_merged = !as.logical(opts_flow$get("split_by_chr"))

  ## if file is available determine whether to split for faster processing
  ## even if the bam files are pre-split its better to explicitily supply the
  ## chr info
  if(split_by_chr){

    ##chrs_info <- get_bam_chrs(tumor_bam)
    chrs_info <- get_fasta_chrs(ref_fasta)
    chrs_prefix <- paste0(bam_prefix, "_", chrs_info) ## bam names
    intervals_opts = paste0(" -L ", chrs_info)             ## interval files

    if(is_merged & length(tumor_bam) > 1){
      stop("multiple bams supplied, expected 1; perhaps is_merged should be FALSE?")

    }else if(!is_merged & length(tumor_bam) == 1){
      stop("single bam supplied, expected multiple (one for each chromosome); perhaps is_merged should be TRUE")
    }
  }else{
  ## dont split
    chrs_prefix = paste0(bam_prefix, ".")
    intervals_opts = ""
  }

  check_args(ignore = "outfile")

  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  mutects <- paste0(chrs_prefix, ".mutect.txt")
  wigs <- paste0(chrs_prefix, ".wig")

  if(!(length(intervals_opts) == length(wigs) | length(intervals_opts) == 1))
    stop("length of intervals should be same as wigs OR 1")

  cmd_mutect <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s --analysis_type MuTect --reference_sequence %s --input_file:tumor %s --input_file:normal %s --out %s  --coverage_file %s %s %s",
                        java_exe, mem_mutect, java_tmp, mutect_jar, ref_fasta, tumor_bam, normal_bam,
                        mutects, wigs,
                        mutect_opts, intervals_opts)
  cmds <- list(mutect = cmd_mutect)

  ## .filter='judgement==KEEP'
  if(split_by_chr){
    merged_mutect = paste0(bam_prefix, "_merged.mutect.tsv")
    cmd_merge = sprintf("flowr ngsflows::merge_sheets x=%s outfile=%s",
                        paste(mutects, collapse = ","), merged_mutect)
    cmds = c(cmds, mutect_merge = cmd_merge)
  }

  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles=list(all = mutects, merged = merged_mutect)))

}


#' Title
#'
#' @param fl a vector of files to be merged
#' @param outfile path to the merged output file
#'
#' @export
merge_sheets <- function(x, outfile, .filter = NA, ...){
  tmp <- lapply(x, function(fl){
    message(".", appendLF = FALSE)
    tab = read_sheet(fl, ...)
    # convert all columns to character, fool proofing
    tab[] <- lapply(x, as.character)

    if(!is.na(.filter)){
      tab2 = dplyr::filter_(tab, .filter)
      return(tab2)
    }else{
      return(tab)
    }
  })

  #mrgd = try(do.call(rbind, tmp))
  # if fails try using dplyr
  mrgd = try(bind_rows(tmp))

  if(!missing(outfile))
    write_sheet(mrgd, outfile)

  invisible(mrgd)
}
