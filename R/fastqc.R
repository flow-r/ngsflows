

#' Fastqc
#'
#' @param fqs 
#' @param fqpath 
#' @param odir 
#' @param fastqc_exe 
#' @param cpu_fastqc 
#' @param casava 
#' @param fastqc_opts 
#'
#' @export
#'
fastqc <- function(fqs, 
         fqpath,
         odir,
         fastqc_exe = get_opts("fastqc_exe"), 
         cpu_fastqc = get_opts("cpu_fastqc"),
         casava = TRUE, 
         fastqc_opts = get_opts("fastqc_opts")
         ){
  
  if(!missing(fqpath))
    fqs = list.files(fqpath, pattern = "fastq.gz", full.names = TRUE, recursive = TRUE)
  
  if(!mean(file.exists(fqs))) 
    stop("Some files do not exist, please check")
  
  ## if input is casava: want to get a summary on fastqs
  if(casava){
    fastqc_opts = c(fastqc_opts, "--casava")
    fqs = paste(fqs, collapse = " ")
  }
  
  fastqc_opts = paste(fastqc_opts, collapse = " ")
  cmds <- sprintf("%s -f fastq -o %s -t %s %s %s",
                  fastqc_exe, odir, cpu_fastqc, fastqc_opts, fqs)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  
  return(list(flowmat = flowmat))
}
attr(fastqc, "type", "module" ) 
