split.names.fastq <- function(files,format="$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"){
  ## process format:
  f <- gsub("\\$samplename\\$","(.*)",format)
  f <- gsub("\\$index\\$","([ATGC]*|NoIndex)",f)
  f <- gsub("\\$lane\\$","([0-9]*)",f)
  f <- gsub("\\$read\\$","([0-9]*)",f)
  f <- gsub("\\$num\\$","([0-9]*)",f)
  ## column names
  out=strsplit(format,"\\$")[[1]]
  cols=out[seq(2,length(out),by=2)]
  repl <- paste("\\",1:length(cols),sep="",collapse=",")
  mat <- gsub(f,repl,basename(files))
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
#' @keywords samplesheet fastq casava
#' @export
#' @examples
#' \dontrun{
#' create_sample_mat(levelipath)
#' }
create_sample_sheet <- function(path, project, subproject, runid, format,
                                pattern = "fastq.gz|fq.gz|fastq|fq", outfile){
  fqs <- unlist(lapply(path, list.files, pattern = pattern,full.names=TRUE,recursive=TRUE))
  if(missing(project)) {project = basename(path); cat("\nDetecting project name:", project)}
  if(missing(subproject)){subproject = substr(project, 1, 2); cat("\nDetecting subproject:", subproject)}
  if(missing(runid)){runid = basename(dirname(dirname(dirname(fqs[1])))) ; cat("\nDetecting runid:", runid)}## runid
  if(missing(outfile)){
    outfile = sprintf("%s_%s_%s_sample_mat.csv", project, subproject, runid);
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
  ## ------- cleanup things
  fq_mat$samplename <- gsub("_","",fq_mat$samplename)
  cat("\nThere are", length(unique(fq_mat$samplename)), "samples in this folder")
  out_basename <- sprintf("%s-%s-%s_%s", project, subproject, fq_mat$samplename, runid)
  sorted_bam <- sprintf("%s_rg.sorted.bam",out_basename)
  recal_bam <- sprintf("%s_rg.sorted.recalibed.bam", out_basename)
  fq_mat <- cbind(fq_mat, out_basename, sorted_bam, recal_bam, runid, project, subproject)
  fq_mat = fq_mat[!grepl("Undetermined", fq_mat$samplename), ] ## remove undermined
  outpath = dirname(outfile)
  if(!file.exists(outpath) & outpath!='.') dir.create(outpath) ## is X exists and not 'blank'
  write.csv(fq_mat, file=outfile, row.names=FALSE)
  return(fq_mat)
}

#' check_fastq_sheet
#' @param mat
#' @export
check_fastq_sheet <- function(mat){
  cat("There are", length(unique(mat$samplename)), "samples in this dataset\n")
  dat_list <- split.data.frame(mat, mat$samplename)
  tmp <- sapply(1:length(dat_list), function(i){
    tmp <- tapply(dat_list[[i]]$files, dat_list[[i]]$read, length)
    if(diff(tmp) != 0) stop("Number of fastq files are not the same in",
                            dat_list[[i]]$samplename[1], "sample")
    return(0)
  })
}

#' read_sample_sheet
#' @param x
#' @export
read_sample_sheet <- function(x){
  ext <- file_ext(x)
  if(ext=="tsv"){
    mat <- read.table(x, as.is=TRUE, sep="\t", header=TRUE)
  }else if(ext=="csv"){
    mat <- read.csv2(x, as.is=TRUE, comment.char = '#', strip.white=TRUE,
                     blank.lines.skip=TRUE, sep=",", header=TRUE)
  }
  check_fastq_sheet(mat)
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
