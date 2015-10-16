
check_params <- function(x){
  needtoexit = FALSE
  for(i in 1:length(x)){
    if(is.null(x[[i]])){
      message("Seems this argument is null, please check: ", names(x[i]))
      needtoexit = TRUE
    }
  }
  if(needtoexit) stop("Some params are not defined properly, please check and resubmit.")
}

#' picard_rg
#' @export
#' @importFrom tools file_path_sans_ext
picard_rg <- function(x, 
                      samplename = get_opts("samplename"),
                      lane = "lane1",
                      ## convert these into get option also, only for this flow
                      seq_platform = "illumina", 
                      center = "MDA",
                      java_exe = get_opts("java_exe"),
                      java_mem = get_opts("java_mem"),
                      java_tmp = get_opts("java_tmp"),
                      picard_dir = get_opts("picard_dir")
){
  
  
  check_args()

  ## make this editable later ....
  rgid = rglb = rgsm = samplename
  rgpu = lane
  
  ## add RG to the orignal bam name
  bamrg_files = sprintf("%s_rg.bam", file_path_sans_ext(x))
  cmds = list(fixrg = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar AddOrReplaceReadGroups INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT",
                              java_exe, java_mem, java_tmp, picard_dir, 
                              x, bamrg_files, rgid, rglb, 
                              seq_platform, rgpu, rgsm, center))
  
  flowmat = to_flowmat(cmds, samplename)
  ret = list(outfiles = bamrg_files, flowmat = flowmat)
  return(ret)
  
}


#' Use picard's MergeSamFiles tool to merge bam/sam files
#' 
#' @description 
#' The resulting file is sorted and index is created for it. 
#' Validation stringency of inputs is kept as lenient.
#' Multi-threading is turned on by default, though in our experience this does
#' not seem to use a lot of threads.
#' 
#' @param x a vectors of files to merge
#' @param mergedbam 
#' @param samplename 
#' @param java_exe 
#' @param java_mem 
#' @param java_tmp 
#' @param picard_dir 
#'
#' @export
picard_merge <- function(x, 
                         mergedbam,
                         samplename = get_opts("samplename"),
                         java_exe = get_opts("java_exe"),
                         java_mem = get_opts("java_mem"),
                         java_tmp = get_opts("java_tmp"),
                         picard_dir = get_opts("picard_dir")){
  

  
  check_args()  
  
  bam_list = paste("INPUT=", x, sep = "", collapse = " ")
  cmds = list(merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar MergeSamFiles %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
                              java_exe, java_mem, java_tmp, picard_dir, bam_list, mergedbam))
  
  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  return(list(outfiles = mergedbam, flowmat = flowmat))
}


#' Title
#'
#' @param x 
#' @param mergedbam 
#' @param samplename 
#' @param java_exe 
#' @param java_mem 
#' @param java_tmp 
#' @param picard_dir 
#'
#' @export
picard_bam_fastq <- function(bam, 
                             samplename = get_opts("samplename"),
                             paired = TRUE,
                             split = FALSE,
                             num_reads = 8000000,
                             split_fq_exe = system.file("scripts/split_fq", package = "ngsflows"),
                             java_exe = get_opts("java_exe"),
                             java_mem = get_opts("java_mem"),
                             java_tmp = get_opts("java_tmp"),
                             picard_dir = get_opts("picard_dir"),
                             bam_fastq_opts = "INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=false RE_REVERSE=true VALIDATION_STRINGENCY=LENIENT"){
  
  check_args()  
  
 
  if(paired){
    fq1 = gsub(".bam$", "_1.fastq", bam)
    fq2 = gsub(".bam$", "_2.fastq", bam)
    fq3 = gsub(".bam$", "_unpaired.fastq", bam)
    fqs = list(fq1 = fq1, fq2 = fq2, fq3 = fq3)
    fqout = sprintf("FASTQ=%s SECOND_END_FASTQ=%s UNPAIRED_FASTQ=%s",
                    fq1, fq2, fq3)
    # we would stream picard into creating splitted fq files. A second flow would read these and start the next step.
    
  }else{
    fq1 = gsub(".bam$", ".fastq", bam)
    fqout = paste0("FASTQ=", fq1)
    fqs = list(fq1 = fq1)
  }

  if(split){
    fq = gsub(".bam$", "", bam)
    fqout = sprintf("FASTQ=/dev/stdout INTERLEAVE=true | bash %s -n %s -f /dev/stdin -o %s",
                    split_fq_exe, as.integer(num_reads), fq)
  }
  
  
  bam_fastq = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar SamToFastq INPUT=%s %s %s",
                                 java_exe, java_mem, java_tmp, picard_dir, bam, bam_fastq_opts, fqout)
  

  cmd = list(bam_fastq = bam_fastq)
  #system(unlist(cmd))
  
  # INPUT is a NAMED list ------
  flowmat = to_flowmat(cmd, samplename)
  return(list(flowmat = flowmat, outfiles = fqs))
}


if(FALSE){
  ## test example
  fqs1 = rep("read1.fq", 20)
  fqs2 = rep("read2.fq", 20)
  out_bwa = bwa(fqs1 = rep("read1.fq", 20),fqs2 = rep("read2.fq", 20))
  
  undebug(to_flowmat)
  #debug(picard_rg);
  out_pic = picard_rg(x = out_bwa$outfiles, samplename = "smp")
  
  
  
  
  
}