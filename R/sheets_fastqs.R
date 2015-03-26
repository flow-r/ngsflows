##
## Author: Sahil Seth
## Date:   2014/09/18
## sseth@mdanderson.org

#' @title split.names.fastq
#' @description split.names.fastq
#' given a format split the files provided. Tried and tested on fastq files
#' @param files
#' @param format
#' @export
split.names.fastq <- function(files,format="$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"){
  ## process format:
  fmt <- gsub("\\$samplename\\$","(.*)",format)
  fmt <- gsub("\\$index\\$","([ATGC]*|NoIndex)", fmt)
  fmt <- gsub("\\$lane\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$read\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$num\\$","([0-9]*)", fmt)
  ## column names
  out = strsplit(format,"\\$")[[1]]
  ## every second one would be a variable
  cols = out[seq(2,length(out), by=2)]
  repl <- paste("\\",1:length(cols),sep="",collapse=",")
  mat <- gsub(fmt, repl,basename(files))
  mat <- data.frame(cbind(do.call(rbind,strsplit(mat,",")),files))
  colnames(mat) <- c(cols,"files")
  return(mat)
}
split_names_fastq=split.names.fastq

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
    if(grepl("_S1.*fastq.gz",fqs[1])){ ## miseq output
      ## ------ casava output
      cat("\nUsing MiSeq naming format")
      format <- "$samplename$_S[0-9]*_L00$lane$_R$read$_$num$.fastq.gz"
    }else if(grepl(".*_([ATGC]*|NoIndex).*L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz",basename(fqs[1]))){
      ## ------ casava output
      cat("\nUsing CASAVA naming format")
      format <- "$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"
    }else{
      stop(c("Looks like we could not understand pattern in names of fastq files",
             print(head(fqs))))
    }
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

#' check_fastq_sheet
#' @param mat
#' @export
check_fastq_sheet <- function(mat, id_column = "sample_id", file_column = "files", read_column = "read"){
  cat("There are", length(unique(mat[, id_column])), "samples in this dataset\n")
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
}

#' read_sample_sheet
#' @param x
#' @export
read_sample_sheet <- function(x, id_column = "sample_id"){
  ext <- tools:::file_ext(x)
  if(ext %in% c("tsv", "txt")){
    mat <- read.table(x, as.is=TRUE, sep="\t", header=TRUE, stringsAsFactors = FALSE,
                      comment.char = '#', strip.white=TRUE, blank.lines.skip=TRUE, quote = "")
  }else if(ext=="csv"){
    mat <- read.csv2(x, as.is=TRUE, sep=",", header=TRUE, stringsAsFactors = FALSE,
                     comment.char = '#', strip.white=TRUE, blank.lines.skip=TRUE, quote = "")
  }
  else if(ext=="xlsx"){
    library(xlsx)
    mat <- xlsx:::read.xlsx2(file = x, sheetName = "sample_sheet", startRow = 2, stringsAsFactors = FALSE)
  }
  else{
    cat("Sorry we do not recognize this file format", ext, "please use tsv, csv or xlsx2 (sheetname: sample_sheet)")
  }
  ### ------ remove blank rows and columns
  mat <- mat[!mat[, id_column] %in% c("", NA), !grepl("^X", colnames(mat))]
  check_fastq_sheet(mat, id_column = id_column)
  return(mat)
}

if(FALSE){

    path = "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-NI56-1"
    runid= "140626_M01692_0046_000000000-A99M8"
    runid="140627_M01692_0047_000000000-A9J0V";
    #debug(create_sample_mat)

    path=sprintf("/scratch/iacs/gcc/leveli/%s",runid)
    fq_mat <- create_sample_mat(path,project="ANDY",subproject="Futreal-AS", runid=runid,
                                outpath="~/flows/ANDY-Futreal-AS")
}
