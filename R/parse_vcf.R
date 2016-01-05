
#' Parse a VCF file into a TSV to be processed and annotated
#'
#' @param x a vcf file to be parsed
#'
#' @export
parse_vcf <- function(x){
  require(data.table)
  require(dplyr)
  require(parallel)

  message("Reading file...")

  hd = scan(x, what = "character", sep = "\n", n = 1000, quiet = TRUE) ##assuming that header is not longer than 100
  hd_row = grep("#CHROM", hd)
  hd = tolower(strsplit(hd[hd_row], "\t")[[1]])
  tab = fread(x, data.table = FALSE, header = FALSE, skip = hd_row, sep = "\t", colClasses = "character")
  colnames(tab) = gsub("#chrom", "chrom", hd)
  tab = tbl_df(tab)

  message("Parsing format columns...")
  samp = colnames(tab)[grep("format", colnames(tab)) + 1]
  cols = splt_vcf_format(x = tab[, samp], format = tab$format, prefix = "")

  message("Parsing info columns...")
  infocols = splt_vcf_info(tab$info)

  colsel = !colnames(tab) %in% c("format", samp, "info")
  message("Assembling data")
  tab = tbl_df(cbind(tab[, colsel], cols, infocols))
  return(tab)
}


splt_vcf_format <- function(x, format, prefix){
  x = as.character(unlist(x))
  format = as.character(unlist(format))
  lst <- lapply(1:length(x), function(i){
    message(".")
    xi = x[i];formati = format[i]
    splt = strsplit(xi, ":")[[1]]
    nms = tolower(strsplit(formati, ":")[[1]])
    names(splt) = paste(prefix, nms, sep = "")
    ret = as.data.frame(t(splt), stringsAsFactors = FALSE)
    return(ret)
  })
  mat = bind_rows(lst)
}


splt_vcf_info <- function(x){
  x = as.character(unlist(x))
  lst <- lapply(1:length(x), function(i){
    #message(i)
    xi = x[i];
    splt = strsplit(xi, ";")[[1]]
    splt = gsub("(.?)\\=(.*)", "\\1\n\n\\2", splt)
    splt = strsplit(splt, "\n\n")
    vals = sapply(splt, tail, 1)
    names(vals) = tolower(sapply(splt, head, 1))
    ret = as.data.frame(t(vals), stringsAsFactors = FALSE)
    return(ret)
  })
  mat = bind_rows(lst)
}


splt_vcf_func <- function(x){
  x = as.character(unlist(x))

  splt_func <- function(i){
    str_replace_all(x[i], pattern = "\\]", "")
    xi = gsub("[{", "", gsub("}]", "", x[i], fixed = TRUE), fixed = TRUE)
    splt = strsplit(xi, ",")[[1]]
    splt = gsub("(.?)\\:(.*)", "\\1\n\n\\2", splt)
    splt = strsplit(splt, "\n\n")
    vals = sapply(splt, tail, 1)
    names(vals) = tolower(sapply(splt, head, 1))
    ret = as.data.frame(t(vals), stringsAsFactors = FALSE)
    return(ret)
  }

  lst <- lapply(1:length(x), function(i){
    y = try(splt_func(i))
    ifelse(class(y)=="try-error", "", y)
  })
  mat = bind_rows(lst)
}



format_vcf_info_ion <- function(x){
  x$freqt = gsub("Freq", "", x$freqt)
  x$freqn = gsub("Freq", "", x$freqn)

}





#' Parse a somatic VCF, with two samples.
#'
#' @param x a vcf file
#' @param samp name of the 'tumor' sample
#' @param ref name of the 'reference' sample
#'
#' @export
#'
parse_somatic_vcf <- function(x, samp, ref){
  message("Reading file...")
  mat = IACSUtil:::readVCF(x)
  mat = tbl_df(as.data.frame(mat,  stringsAsFactors = FALSE))
  colnames(mat) = tolower(colnames(mat))
  ## assume first column
  samp = colnames(mat)[grep("format", colnames(mat)) + 1]
  ref = colnames(mat)[grep("format", colnames(mat)) + 2]
  message("Parsing format columns...")
  tcols = splt_vcf_format(x = mat[,samp], format = mat$format, prefix = "t_")
  ncols = splt_vcf_format(x = mat[,ref], format = mat$format, prefix = "n_")
  message("Parsing info columns...")
  infocols = splt_vcf_info(mat$info)
  message("Parsing func...")
  #funccols = splt_vcf_func(mat$func)
  #infocols[1:100,] %>% View()
  mat$t_sample = samp
  mat$n_sample = ref
  colsel = !colnames(mat) %in% c("format", samp, ref, "info")
  message("Assembling data")
  #mat = tbl_df(cbind(mat[, colsel], tcols, ncols, infocols, funccols))
  mat = tbl_df(cbind(mat[, colsel], tcols, ncols, infocols))
  return(mat)
}

if(FALSE){

  source('~/Dropbox/public/github_ngsflows/R/parse_vcfs.R')
  x='/rsrch2/iacs/iacs_dep/sseth/flowr/runs/flowr_test/fastq_haplotyper-MS132-20150824-16-37-58-XScJT0OZ/tmp/GLizee-Pancreatic-MS132-MP013normalDNA.recalibed_1.haplotyper.vcf'
  x2 = parse_vcf(x)

}
