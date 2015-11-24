##
## Author: Sahil Seth
## Date:   2014/09/18
## Edited: 2015/07/15
## major edits
## sseth@mdanderson.org

detect_fq_format2 <- function(x){
  if(grepl("_S1.*fastq.gz", x[1])){ ## miseq output
    ## ------ casava output
    message("Using MiSeq/bcl2fastq (HiSeq3000) naming format")
    format <- "{{samplename}}_S[0-9]*_L00{{lane}}_R{{read}}_{{num}}.fastq.gz"
  }else if(grepl(".*_([ATGC]*|NoIndex).*L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz",basename(x[1]))){
    ## ------ casava output
    message("Using CASAVA 1.8 naming format")
    format <- "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz"
  }else{
    stop(c("Looks like we could not understand pattern in names of fastq files",
           print(head(x))))
  }
}


#' split_names_fastq2
#' @importFrom flowr whisker_render
split_names_fastq2 <- function(x, format = "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz"){

  ## --- regex pattern for each piece
  lst_patterns = list(
    samplename = "(.*)",
    index = "([ATGC]*|NoIndex)",
    lane = "([0-9]*)",
    read = "([0-9]*)", ## ideally would be 1 OR 2
    num = "([0-9]*)")

  ## --- replace and get final pattern !
  tmp = whisker_render(format, lst_patterns)
  fmt = tmp$out
  vars = tmp$vars

  ## get the parse matrix
  repl <- paste("\\",1:length(vars), sep="",collapse=",")
  mat <- gsub(fmt, repl, basename(x))

  ## longer, but more robust
  mat = lapply(1:length(mat), function(i){
    out = strsplit(mat[i], ",")[[1]]
    add_vars = length(vars) - length(out)
    if(add_vars > 0){
      out = c(out, rep(NA, add_vars))
      warning("there was a issue parsing this filename: ", basename(x[i]))
    }
    return(out)
  })
  warnings()

  mat <- do.call(rbind, mat)
  mat <- data.frame(mat, x, stringsAsFactors = FALSE)

  ## if all files were not parsed properly
  ## show extra message, to help debug
  ## This does not catch, when some file are parsed improperly.
  if(ncol(mat) != length(c(vars, "file"))){
    message("Final evaluated regular expression: ", fmt, " ---> ", repl)
    message("Format was not able to split fastq names properly.")
    print(head(mat))
    stop("Exiting create_fq_sheet...")
  }

  ## if all is well
  colnames(mat) <- c(vars, "file")
  return(mat)
}

#' check_fastq_sheet
#'
#' Checks the nuber of samples, and that number of files in each sample should be
#' same for each read (if paired end)
#' @param x a data.frame with details regarding fastq files
#'
#' @export
check_fastq_sheet <- function(x,
                              id_column = "samplename",
                              file_column = "file",
                              read_column = "read"){
  mat = x
  message("There are ", length(unique(mat[, id_column])), " samples in this dataset")
  dat_list <- split.data.frame(mat, mat[, id_column])
  if(length(mat$files) > 0){
    tmp <- sapply(1:length(dat_list), function(i){
      tmp <- tapply(dat_list[[i]][, file_column], dat_list[[i]][, read_column], length)
      ## check in case of paired end
      if(length(tmp) > 1){
        if(diff(tmp) != 0) stop("Number of fastq files are not the same in",
                                dat_list[[i]][, id_column][1], "sample")
      }
      return(0)
    })}
  return(x)
}




options(
  ngs_fq_ext = "fastq.gz"
)

#' create a sheet of fastq
#' @param x path to a fastq folder
#' @param ext extensions to look for. A regex to search for fastq files
#' @param format auto detect. Pattern specify pattern acceptable to split_fastq_names. If missing will detect hiseq and miseq
#' @param sample_prefix A prefix to add to all sample names, run, project etc.....
#' @export
create_fq_sheet <- function(x,
                            ext = opts_flow$get("fastq_extension"),
                            format = opts_flow$get("fastq_format"),
                            sample_prefix = ""){

  ## get the extension to use
  if(is.null(ext)) ext = "fastq.gz|fastq"
  if(!grepl("$", ext, fixed = TRUE))
    ext = paste0(ext, "$") ## add a dollar

  ## fetch the files to parse
  fqs <- unlist(lapply(x, list.files, pattern = ext, full.names=TRUE,recursive=TRUE))
  if(length(fqs) == 0)
    stop("no fastq files found. Please check the folder provided.")

  ## get the format to use
  if(missing(format))
    format <- detect_fq_format2(fqs)
  if(is.null(format))
    format <- detect_fq_format2(fqs)

  fqmat <- split_names_fastq2(fqs, format)

  ## if one needs to change the samplename to include anything in this
  fqmat$samplename = paste0(sample_prefix, fqmat$samplename)

  return(fqmat)
}



if(FALSE){
  x = "/rsrch1/genetics/htep/hiseq/150707_ST-J00106_0022_AH35LGBBXX/Unaligned/YLD-NS-KO_RNA109"
  fqmat = create_fq_sheet(x)
  fqmat = check_fastq_sheet(fqmat)

}

########################## ------------------- functions below are not used
#      are here for archival purposes only





if(FALSE){

  path = "/rsrch1/genetics/htep/hiseq/150707_ST-J00106_0022_AH35LGBBXX/Unaligned/YLD-NS-KO_RNA109"
  runid= "140626_M01692_0046_000000000-A99M8"
  runid="140627_M01692_0047_000000000-A9J0V";
  #debug(create_sample_mat)

  path=sprintf("/scratch/iacs/gcc/leveli/%s",runid)
  fq_mat <- create_sample_mat(path,project="ANDY",subproject="Futreal-AS", runid=runid,
                              outpath="~/flows/ANDY-Futreal-AS")
}





detect_fq_format <- function(x){
  .Deprecated("detect_fq_format2")
  if(grepl("_S1.*fastq.gz", x[1])){ ## miseq output
    ## ------ casava output
    message("\nUsing MiSeq naming format")
    format <- "$samplename$_S[0-9]*_L00$lane$_R$read$_$num$.fastq.gz"
  }else if(grepl(".*_([ATGC]*|NoIndex).*L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz",basename(x[1]))){
    ## ------ casava output
    message("\nUsing CASAVA naming format")
    format <- "$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"
  }else{
    stop(c("Looks like we could not understand pattern in names of fastq files",
           print(head(x))))
  }
}



#' Creates a sample sheet file names in the provided folder
#'
#' This function would check for files ending in fastq.gz, fq.gz, fastq, fq.
#'
#' @param path path to the fastq files
#' @param project name of the project.  \emph{optional}
#' @param subproject name of the subproject \emph{optional}
#' @param runid name of the flowcell this data is from \emph{optional}
#' @param outfile name of the output csv files \emph{optional}
#' @param format the format for names of fastq files, we have defaults for CASAVA and miSeq
#' @param pattern extensions this function will look for \emph{optional}
#' @param fix.names change the sample names such that they are acceptable as column names for R
#' @keywords samplesheet fastq casava
#'
#' @export
#' @examples
#' \dontrun{
#' create_sample_mat(levelipath)
#' }
create_sample_sheet <- function(path, project, subproject, runid, format,
                                fix.names = FALSE,  fix.names.char = "-",
                                out_sep = c("\t", ","),
                                include.undetermined = FALSE,
                                pattern = "fastq.gz|fq.gz|fastq|fq", outfile){

  message("Fetching path(s) for fastq files...\n")
  .Deprecated("create_fq_sheet")
  fqs <- unlist(lapply(path, list.files, pattern = pattern,full.names=TRUE,recursive=TRUE))
  if(!include.undetermined) fqs <- grep("Undetermined", fqs, value = TRUE, invert = TRUE)
  if(length(fqs) == 0) stop("No fastq files detected in this folder\n")
  if(missing(project)) {project = basename(path); cat("\nDetecting project name:", project)}
  if(missing(subproject)){subproject = substr(project, 1, 2); cat("\nDetecting subproject:", subproject)}
  if(missing(runid)){runid = basename(dirname(dirname(dirname(fqs[1])))) ; cat("\nDetecting runid:", runid)}## runid
  out_sep = match.arg(out_sep)
  if(missing(outfile)){
    outfile = sprintf("%s_%s_%s_sample_mat.%s", project, subproject, runid,
                      switch(out_sep,
                             "," = "csv",
                             "\t" = "tsv"))
    cat("\nDetecting outfile:", outfile)
  }## folder for samplemat
  if(missing(format)){
    format <- detect_fq_format(fqs[1])
  }
  fq_mat <- split_names_fastq(files = fqs, format = format)
  if(fix.names){
    fq_mat$sample_id_orig = fq_mat$sample_id
    fq_mat$sample_id = fix_names(fq_mat$sample_id, char = fix.names.char)
    if(fix.names.char == ".") ## . opt style 2, good for data.frames
      fq_mat$sample_id = make.names(fq_mat$sample_id)
    ## ------- cleanup things: use _ to seperate out other things
  }
  cat("\nThere are", length(unique(fq_mat$samplename)), "samples in this folder")
  out_basename <- sprintf("%s-%s-%s_%s", project, subproject, fq_mat$samplename, runid)
  ## sorted_bam <- sprintf("%s_rg.sorted.bam",out_basename)
  ## recal_bam <- sprintf("%s_rg.sorted.recalibed.bam", out_basename)
  fq_mat <- cbind(fq_mat, out_basename, runid, project, subproject)
  fq_mat = fq_mat[!grepl("Undetermined", fq_mat$samplename), ] ## remove undermined
  outpath = dirname(outfile)
  if(!file.exists(outpath) & outpath!='.') dir.create(outpath) ## is X exists and not 'blank'
  write.table(fq_mat, file = outfile, row.names=FALSE, sep = out_sep, quote = FALSE)
  return(fq_mat = fq_mat)
}


#' @title split.names.fastq
#' @description Given a format split the files provided. Tried and tested on fastq files
#'
#' @param files a character vector of filenames, to be parsed
#' @param format the regex type format to be used to extract information from filenames. See default format as an example.
#'
#' @export
split_names_fastq <- function(files,format="$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"){
  ## process format:
  .Deprecated("split_names_fastq2")
  fmt <- gsub("\\$samplename\\$","(.*)",format)
  fmt <- gsub("\\$index\\$","([ATGC]*|NoIndex)", fmt)
  fmt <- gsub("\\$lane\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$read\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$num\\$","([0-9]*)", fmt)
  ## -- column names
  out = strsplit(format,"\\$")[[1]]
  ## every second one would be a variable
  cols = out[seq(2,length(out), by=2)]
  repl <- paste("\\",1:length(cols),sep="",collapse=",")
  mat <- gsub(fmt, repl,basename(files))
  mat <- data.frame(cbind(do.call(rbind,strsplit(mat,",")),files))
  colnames(mat) <- c(cols,"files")
  return(mat)
}
split.names.fastq=split_names_fastq
