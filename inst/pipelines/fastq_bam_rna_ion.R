## RNA-Seq pipeline
## http://ioncommunity.lifetechnologies.com/docs/DOC-7062
## remove adaptor sequences
## cutadapt -m 16 -b GGCCAAGGCG -o adaptorTrim.fastq input.fastq
## align with tophat2
## tophat2 -p 12 --keep-fasta-order --GTF known_genes.gtf \
## hg19_bowtie2_index adaptorTrim.fastq
## convert unmapped to fastq from bam
## bam2fastq -o unmapped.fastq -q unmapped.bam
## align the unmapped with bowtie2
## bowtie2 --local --very-sensitive-local -p 8 --mm -x hg19_bowtie2_index \
## -U unmapped.fastq  | samtools view -uhS -F4 - | samtools sort - unmapped_remap
## merge bam files with Picard module MergeSamFiles
## java -jar MergeSamFiles USE_THREADING=true MSD=true AS=true I=accepted_hits.bam I=unmapped_remap.bam O=aligned.bam


## http://blog.sbgenomics.com/ion-proton-rna-seq-alignment/


set_opts(
  cutadapt_exe = "module load cutadapt;cutadapt",
  tophat2_exe = "module load tophat2;tophat2",
  cutadapt_exe = "module load cutadapt;cutadapt",
  bam2fastq_exe = "bam2fastq"
)



aln_rna_ion <- function(fq,
                        samplename = "mysamp",
                        cutadapt_exe = get_opts("cutadapt_exe"),
                        adap_seq = "GGCCAAGGCG",
                        tophat2_exe = get_opts("tophat2_exe"),
                        star_exe = get_opts("star"),
                        bam2fastq_exe = get_opts("bam2fastq_exe"),
                        bowtie2_exe = get_opts("bowtie2_exe"),
                        java_exe = get_opts("java_exe"),
                        picard_dir = get_opts("picard_dir"),
                        gtf = get_opts("ref_gtf"),
                        ref_bowtie2 = get_opts("ref_bowtie2"),
                        ref_star = get_opts("ref_bowtie2")
                        
                        ){
  
  assert_that(is.character(fq))
  assert_that(is.character(cutadapt_exe))
  assert_that(is.character(adap_seq))
  assert_that(is.character(tophat2_exe))
  assert_that(is.character(bam2fastq_exe))
  assert_that(is.character(bowtie2_exe))
  assert_that(is.character(java_exe))
  assert_that(is.character(picard_dir))
  assert_that(is.character(gtf))
  assert_that(is.character(ref_bowtie2))
  
  ## assert that none of the arguments are null
  check_args()
  
  fq_adap = gsub("fastq$", "adapt_trim.fastq", fq)
  trim = sprintf("%s -m 16 -b %s -o %s %s",
                 cutadapt_exe, adap_seq, fq_adap, fq)
  
  tophat = sprintf("tophat2 -p 12 --keep-fasta-order --GTF %s %s %s", 
                   gtf, ref_bowtie2, fq_adap)
  
  #star = star(fq1 = fq)

  ## convert unmapped to fastq from bam
  bam2fastq = sprintf('bam2fastq -o unmapped.fastq -q unmapped.bam')
  
  ## align the unmapped with bowtie2
  bowtie = sprintf("bowtie2 --local --very-sensitive-local -p 8 --mm -x %s -U unmapped.fastq | samtools view -uhS -F4 - | samtools sort - unmapped_remap",
                   bowtie2_exe, ref_bowtie2)
  
  ## merge bam files with Picard module MergeSamFiles
  outbam = gsub("fastq$", "bam", fq)
  merge = sprintf("%s -jar %s/MergeSamFiles.jar USE_THREADING=true MSD=true AS=true I=accepted_hits.bam I=unmapped_remap.bam O=%s",
                  java_exe, picard_dir, outbam)
  
  cmds = list(trim = trim, tophat = tophat, bam2fastq = bam2fastq, bowtie = bowtie, merge = merge)
  
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat = flowmat, outfiles = outbam))
}


if(FALSE){
  
  require(assertthat)
  devtools::load_all()
  

  out = load_opts("inst/conf/ngsflows_mda.conf")
  set_opts(star_opts = "", star_cpu = 8)
  fq = "/rsrch1/iacs/ngs_runs/Sarcomatoid/rna_seq/fastq/downloads/RNA_Barcode_None_001.R_2014_09_16_13_13_35_user_PRO-99-Sarcomatoid_RNA_Seq_2738-T.fastq"
  #debug(aln_rna_ion)
  out = aln_rna_ion(fq)
  
  fq2 = "/rsrch1/iacs/ngs_runs/Sarcomatoid/rna_seq/fastq/downloads/RNA_Barcode_None_001.R_2014_09_16_13_13_35_user_PRO-99-Sarcomatoid_RNA_Seq_2738-T.fastq"
  out = star(fq1 = fq2)
  cat(paste(out$flowmat$cmd, collapse = "\n"))
  
  
  bam = "/rsrch2/iacs/ngs_runs/1411_sarcomatoid/bams-2015.07.24/RNA_Barcode_None_user_PRO-100-Sarcomatoid_RNA_Seq_2738-N_230.bam"
  
}

aln_rna_bam_ion <- function(bam){
  
  #debug(star)
  out = star(bam = bam)
  
}


