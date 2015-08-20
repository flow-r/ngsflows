## cpu_markdup=4, java_mem_markdup= "-Xmx8g",
## cpu_target=16,  java_mem_target= "-Xmx32g",
## cpu_realign=4,  java_mem_realign= "-Xmx4g", ## scatter 8 per node
## cpu_baserecalib=4,  java_mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
## cpu_printreads=4,  java_mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
#' Flow following Broad's best practices for variant calling, starting from sorted bam

#' @title Flow following Broad's best practices for variant calling, starting from sorted bam
#' @description This function provides a list of commands to run.
#' CPU's refer to number of cores reserved the job. We have a seperate param to control the number of threads the command would use. Example. On a node with 10GB and 10 cores. We may only reserve 10, but use only 2 threads for a 5GB memory command.
#' @param inbam
#' @param outbam
#' @param q_obj is provided output is a flow object, else a list of commands to run
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' @param cpu_markdup cpus to be reserved, default 4
#' @param java_mem_markdup memory provided to java
#' @param cpu_target
#' @param java_mem_target
#' @param cpu_realign
#' @param java_mem_realign
#' @param cpu_baserecalib
#' @param java_mem_baserecalib
#' @param cpu_printreads
#' @param java_mem_printreads
#' @param gatk_jar
#' @param picard_dir
#' @param reffa
#' @param gatk_target_opt
#' @param gatk_realign_opt
#' @param gatk_baserecalib_opt
#' @param printreads_opt
#' @export
bam_preprocess <- function(inbam, outbam, q_obj,
                           java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                           java_tmp = "tmp",
                           cpu_markdup = 4, java_mem_markdup= "-Xmx8g",
                           cpu_target = 16,  java_mem_target= "-Xmx32g",
                           nt_target = 16, ## reserve 16 use 16
                           cpu_realign = 8,  java_mem_realign= "-Xmx4g", ## scatter 8 per node
                           cpu_baserecalib = 4,  java_mem_baserecalib= "-Xmx4g", ## scatter 8 per node nct=8
                           cpu_printreads = 4,  java_mem_printreads= "-Xmx4g", ## scatter 8 per node nct=8
                           gatk_jar = "/scratch/rists/hpcapps/x86_64/gatk/3.1-1/GenomeAnalysisTK.jar",
                           picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                           reffa = "/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta",
                           gatk_target_opt = "-known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/mills_and_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf -known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores",
                           gatk_realign_opt = "-known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/mills_and_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf -known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores",
                           gatk_baserecalib_opt  = "-knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/dbsnp/dbsnp_138.b37.excluding_sites_after_129.vcf -knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/mills_and_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores -nct 8",
                           printreads_opt = "-allowPotentiallyMisencodedQuals -nct 8"){
  chrs_info <- get_bam_chrs(bam)
  bam_prefix <- gsub(".bam", "", basename(inbam))
  chrs_prefix <- paste0(gsub(".bam", "_", basename(inbam)), chrs_info[, 1])
  
  dedupbam <- gsub(".bam",".marked.bam",basename(inbam))
  metricsfile <- gsub(".bam",".marked.metrics",basename(inbam))
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/MarkDuplicates.jar INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORTED=true   VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                         java_exe, java_mem_markdup, java_tmp,picard_dir, 
                         inbam, dedupbam, metricsfile)
  ## ------------ realign
  intervalsfile <- paste0(bam_prefix, ".realign.intervals")
  realignedbams <- paste0(chrs_prefix ,".realigned.bam")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s",
                        java_exe, java_mem_target, java_tmp, gatk_jar, reffa, dedupbam, 
                        intervalsfile, nt_target, gatk_target_opt)
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s",
                         java_exe, java_mem_realign, java_tmp, gatk_jar, reffa, inbam, 
                         intervalsfile, dedupbam, gatk_realign_opt)
  ## ------------ base recalibration
  recalibbam <- gsub(".bam",".recalibed.bam",basename(inbam))
  recalibtabfile <- gsub(".bam",".recalib.tab",basename(inbam))
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s %s",
                             java_exe,java_mem_baserecalib,java_tmp, gatk_jar,reffa,realignedbam,recalibtabfile,
                             gatk_baserecalib_opt)
  cmd_printreads <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s %s",
                            java_exe,java_mem_printreads, java_tmp, gatk_jar,reffa,realignedbam,recalibtabfile,recalibbam,
                            printreads_opt)
  cmds <- c(cmd_markdup=cmd_markdup, cmd_target=cmd_target, cmd_realign=cmd_realign,
            cmd_baserecalib=cmd_baserecalib, cmd_printreads=cmd_printreads)
  jobs <- NA; fobj <- NA
  ## if queue object is not given output commands
  if(!missing(q_obj)){
    flowname="preprocess"
    j_obj_markdup <- job(q_obj=q_obj, name="markdup", cmds=cmd_markdup,  cpu=cpu_markdup,
                         submission_type="serial")
    j_obj_target <- job(q_obj=q_obj, name="target", cmds=cmd_target,  cpu=cpu_target,
                        previous_job="markdup", dependency_type="serial")
    j_obj_realign <- job(q_obj=q_obj, name="realign", cmds=cmd_realign,  cpu=cpu_realign,
                         previous_job="target", dependency_type="serial")
    j_obj_baserecalib <- job(q_obj=q_obj, name="baserecalib", cmds=cmd_baserecalib,  cpu=cpu_baserecalib,
                             previous_job="realign", dependency_type="serial")
    j_obj_printreads <- job(q_obj=q_obj, name="printreads", cmds=cmd_printreads,  cpu=cpu_printreads,
                            previous_job="baserecalib", dependency_type="serial")
    jobs <- list(j_obj_markdup, j_obj_target, j_obj_realign, j_obj_baserecalib,j_obj_printreads)
    flow_desc = gsub("$.bam","",basename(inbam))
    fobj <- flow(jobs=jobs,name=flowname, mode="scheduler", desc=flow_desc)
  }
  return(list(cmds=cmds, flow=fobj))
}

get_bam_chrs <- function(bam, samtools_exe = "samtools") {
  cmd <- sprintf("%s view -H %s | grep  '\\@SQ' | cut -f 2,3",
                 samtools_exe, bam)
  chrs_info <- system(cmd, intern = TRUE)
  chrs_info <- do.call(rbind, lapply(chrs_info, function(x){
    y = strsplit(x, "\t|:")[[1]]
    y[c(2,4)]
  }))
  return (chrs_info)
}


if(FALSE){
  bam="/scratch/iacs/gcc/levelii/140829_SN208_0523_BC5FKTACXX/TCGA-VD-AA8P-10A-01D-A40E_140829_SN208_0523_BC5FKTACXX_s_8_rg.sorted.bam"
  ##testing on
  ## /scratch/iacs/iacs_dep/sseth/flows/AdenoidCysticSanger/PD3176a/merge-preproc-0dd3b13a-aef5-4cdc-bee6-15db06f2dc71/tmp
}
