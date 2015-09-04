
#' Title
#'
#' @param x 
#' @param y 
#' @param samplename 
#' @param outfile 
#' @param is_merged 
#' @param split_by_chr 
#'
#' @export
#'
#' @examples
#' x = "tumor.bam"
#' y = "normal.bam"
#' 
#' out = mutect(x, y, is_merged = TRUE)
#' 
#' x = "tumor.bam"
mutect <- function(x, y,
                   samplename = get_opts("samplename"),
                   outfile,
                   is_merged = TRUE,
                   split_by_chr = TRUE, 
                   
                   java_exe = get_opts("java_exe"),
                   java_tmp = get_opts("java_tmp"),

                   mutect_jar = get_opts('mutect_jar'),
                   
                   cpu_mutect = 1,
                   mem_mutect = "-Xmx8g",
                   
                   ref_fasta = get_opts('ref_fasta'),
                   
                   mutect_opts = get_opts('mutect_opts')

                   ){

    ## determine output file name
    if(missing(outfile))
        bam_prefix <- gsub(".bam", "", basename(x))
    else
        bam_prefix <- gsub(".bam", "", basename(outfile))

    ## if file is available determine whether to split for faster processing
    if(split_by_chr & is_merged){
        ##chrs_info <- get_bam_chrs(x)
        chrs_info <- get_fasta_chrs(ref_fasta)
        chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_") ## bam names
        intervals_opts = paste0(" -L ", chrs_info)             ## interval files

    }else if(split_by_chr & !is_merged){
        chrs_prefix = bam_prefix
        intervals_opts = ""

    }else{
        chrs_prefix = bam_prefix
        intervals_opts = ""
    }

    check_args(ignore = "outfile")

    pipename = match.call()[[1]]
    message("Generating a ", pipename, " flowmat for sample: ", samplename)

    mutects <- paste0(chrs_prefix, ".mutect.txt")
    wigs <- paste0(chrs_prefix, ".wig.txt")

    cmds_mutect <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s --analysis_type MuTect --reference_sequence %s --input_file:tumor --input_file:normal --out %s  --coverage_file %s %s %s",
                           java_exe, mem_mutect, java_tmp, mutect_jar, ref_fasta, x, y,
                           mutects,
                           mutect_opts, inverval_opts)

    cmds <- list(mutect = cmd_mutect)

    flowmat = to_flowmat(cmds, samplename = samplename)
    return(list(flowmat=flowmat, outfile=recalibbams))

}
