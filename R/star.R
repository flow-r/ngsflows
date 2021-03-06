## in-house



star.2.4.0j <- function(fq1, fq2,
  samplename = "samp1",
  out_path, out_prefix,
  star_exe = opts_flow$get("star_exe"),
  ref_fasta = opts_flow$get("ref_fasta"),
  ref_star = opts_flow$get("ref_star"),
  star_opts = opts_flow$get("star_opts"),
  star_cpu = opts_flow$get("star_cpu")){

  .Deprecated("2.4.2a")

  ## detect gz and change options
  if(grepl("\\.gz$", fq1)){
    star_opts = paste0(star_opts, " --readFilesCommand zcat")
  }

  if(missing(fq2))
    fq = fq1
  else
    fq = c(fq1, fq2)

  chkfq <- chk_fq(fqs1 = fq1, fqs2 = fq2)

  if(missing(out_prefix))
    out_prefix = gsub(chkfq$ext, "", basename(fq1))

  prefix1 <- paste0(out_prefix,"initial",sep='') ## prefix
  prefix2 <- paste0(out_prefix,"star", sep='') ## prefix
  ref2_star <- paste0(out_prefix, 'genomedir')

  ## -- files generated by star
  sjtab <- paste0(prefix1, "_SJ.out.tab")
  outsam <- paste0(prefix1, "Aligned.out.sam")
  sjtab2 <- paste0(prefix2, "_SJ.out.tab")
  outsam2 <- paste0(prefix2, "Aligned.out.sam")

  cmd_align1 <- sprintf("%s --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s_ %s",
    star_exe, ref_star, fq, star_cpu, prefix1, star_opts)

  cmd_index = sprintf("mkdir %s; %s --runMode genomeGenerate --genomeFastaFiles %s --sjdbFileChrStartEnd %s --genomeDir %s %s",
    ref2_star, star_exe, ref_fasta, sjtab, ref2_star, star_opts)

  cmd_align2 <- sprintf("%s --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s_ %s",
    star_exe, ref2_star, fq, star_cpu, prefix2, star_opts)

  cmds = list(align1 = cmd_align1, index = cmd_index, align2 = cmd_align2)

  flowmat = to_flowmat(cmds, samplename)
  return(list(flowmat = flowmat, outfiles = c(sjtab2, outsam2)) )
}





## input from ion-torrent is preprocessed:
## @PG     ID:bc   PN:BaseCaller   VN:4.2-14/88431 CL:BaseCaller --barcode-filter 0.01 --barcode-filter-minreads 10 --keypass-filter on --phasing-residual-filter=2.0 --num-unfiltered 1000 --barcode-filter-postpone 1 --barcode-mode 1 --barcode-cutoff 2 --input-dir=sigproc_results --librarykey=TCAG --tfkey=ATCG --run-id=V4IOG --output-dir=basecaller_results --block-col-offset 0 --block-row-offset 7992 --datasets=basecaller_results/datasets_pipeline.json --trim-adapter ATCACCGACTGCCCATAGAGAGGCTGAGAC

#' @rdname star
#' @usage star.2.4.2a
star.2.4.2a <- function(fq1, fq2, bam,
  samplename = "samp1",
  out_path,
  out_prefix,
  ref_fasta = opts_flow$get("ref_fasta"),
  ref_star = opts_flow$get("ref_star"),
  star_exe = opts_flow$get("star_exe"),
  star_opts = opts_flow$get("star_opts"),
  star_cpu = opts_flow$get("star_cpu")
  ){


  assert_that(is.character(ref_fasta))
  assert_that(is.character(ref_star))
  assert_that(is.character(star_exe))
  assert_that(is.character(star_opts))
  assert_that(!is.null(star_cpu)) ## not numeric, might be {{{CPU}}}


  if(!missing(fq1)){
    if(!missing(bam))
      stop("Supply either fq OR bam as input, not both. Refer to STAR's manual for more details.")

    ## detect gz and change options
    if(grepl("\\.gz$", fq1)){
      star_opts = paste0(star_opts, " --readFilesCommand zcat")
    }

    if(missing(fq2))
      fq = fq1
    else
      fq = c(fq1, fq2)

    chkfq <- chk_fq(fqs1 = fq1, fqs2 = fq2)

    if(missing(out_prefix)) ## guess out_prefix
      out_prefix = gsub(chkfq$ext, "", basename(fq1))

    cmd_star <- sprintf("%s  --twopassMode Basic --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s_ %s",
      star_exe, ref_star, fq, star_cpu, out_prefix, star_opts)
  }


  if(!missing(bam)){

    if(missing(out_prefix)) ## guess out_prefix
      out_prefix = gsub(".bam", "", basename(bam))

    ## notes:
    ## EXITING because of fatal INPUT ERROR: at the moment --runMode inputFromBAM only works with --outWigType bedGraph OR --bamRemoveDuplicatesType Identical
    cmd_star <- sprintf("%s  --twopassMode Basic --genomeDir %s --inputBAMfile %s --runMode inputAlignmentsFromBAM --outWigType bedGraph --runThreadN %s --outFileNamePrefix %s_ %s",
      star_exe, ref_star, bam, star_cpu, out_prefix, star_opts)
  }

  if(missing(bam) & missing(fq1))
    stop("Both inputs, bam and fastq are missing !")

  flowmat = to_flowmat(list(star = cmd_star), samplename)
  return(list(flowmat = flowmat, outfiles = c(out_prefix)) )


}



## STAR DOC

#' @title A flowr wrapper for STAR aligner
#'
#' @description For more details refer to \href{https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf}{STAR's manual}.
#' Starting 2.4.2a, star supports a twoPass mode, thus simplifying this pipeline.
#'
#' @usage
#' star(fq1, fq2, bam,
#'  samplename = "samp1", out_path, out_prefix,
#'  ref_fasta = opts_flow$get("ref_fasta"), ref_star = opts_flow$get("ref_star"),
#'  star_exe = opts_flow$get("star_exe"), star_opts = opts_flow$get("star_opts"),
#'  star_cpu = opts_flow$get("star_cpu"))
#'
#' @param fq1 single fastq file
#' @param fq2 second read fastq file
#' @param bam alternatively use of bam as input, instead of fqs
#' @param samplename name of the sample
#' @param out_path output path
#' @param out_prefix output prefix
#' @param star_exe path to star executable
#' @param ref_fasta path to reference genome in fasta format
#' @param ref_star path to star's reference index
#' @param star_opts additional options passed onto star. This is a simple character vector appended to the commandline.
#' @param star_cpu a integer specifying the number of cores to be used
#' @export
star <- star.2.4.2a



opts_flow$set(
  star_opts = "--sjdbOverhang 75",
  star_cpu = "{{{CPU}}}"
)

## would get from flow_def

#' star_index
#' @description a helper function to create a index, incase it does not exits
star_index <- function(
  ref_star = opts_flow$get("ref_star"),
  star_exe = opts_flow$get("star_exe"),
  ref_fasta = opts_flow$get("ref_fasta"),
  star_cpu = opts_flow$get("star_cpu"),
  ref_gtf = opts_flow$get("ref_gtf")){

  sprintf("%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %s",
    star_exe, ref_star,  ref_fasta, star_cpu)

}


if(FALSE){
  ## get a star index


  #debug(star_pipe)
  ref_fasta = '/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta'
  ref_star = '/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/STAR/2.3.0e'

  fq1 = "/rsrch2/iacs/ngs_runs/1411_sarcomatoid/plugin_out/downloads/1711-N.RNA_Barcode_None_001.user_PRO-125-Sarcomatoid_RNA_Seq_1711-N.fastq"
  out = star(fq1, star_cpu = 8, ref_fasta = ref_fasta, ref_star = ref_star)


}


## example star pipeline for ion-torrent
if(FALSE){

  undebug(whisker:::parseTemplate)
  debug(params:::parse_conf)

  load_opts(fetch_conf("ngsflows_mda.conf"))
  opts_flow$get()

  cmdindex = star_index(ref_gtf = opts_flow$get("ref_hg19_gtf"))


  bam = "/rsrch2/iacs/ngs_runs/1411_sarcomatoid/bams-2015.07.24/RNA_Barcode_None_user_PRO-64-Sarcomatoid_RNA_Seq_768-T_190.bam"
  out = star(bam, ref_fasta = ref_fasta, ref_star = ref_star, star_opts = "--sjdbOverhang 100")
  out

}
