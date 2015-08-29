
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
#     if(missing(samplename))
#       stop("this function needs a samplename.")

    ## -----  check if any of the params are null
#     params = lapply(names(formals()), function(zzz) get(zzz))
#     names(params) = names(formals())
#     check_params(params)
    
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
  
  ## -----  check if any of the params are null
#   params = lapply(names(formals()), function(zzz) get(zzz))
#   names(params) = names(formals())
#   check_params(params)

  check_args()  
  
  bam_list = paste("INPUT=", x, sep = "", collapse = " ")
  cmds = list(merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar MergeSamFiles %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
    java_exe, java_mem, java_tmp, picard_dir, bam_list, mergedbam))

  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  return(list(outfiles = mergedbam, flowmat = flowmat))
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