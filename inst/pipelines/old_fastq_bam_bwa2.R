## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------


## ----------------------- d e p r e c i a t e d ---------------------------------------------- ##

as.c=as.character


flow_aln_merge <- function(fqs1, fqs2, out_bam, sample="sample",
                           out_bampath="bams", se.or.pe,
                           runid="runid", project="project", subproject="subproject", lane="0", ## RG tags
                           platform  = 'illumina',center = 'IACS', ## RG tags
                           bwa_exe = "/scratch/iacs/tmp/bwa-0.7.9a/bwa",
                           samtools_exe = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools",
                           picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                           java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",java_mem="-Xmx4g",
                           reflib = "/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/bwa/0.7.5a/Homo_sapiens_assembly19.fasta",
                           bwa_aln_opts = "-l 50 -n 6 -k 4 -t 6",
                           bwa_sampe_opts="",
                           q_type="torque",q_queue="iacs",
                           cpu_aln=12,cpu_sampe=1,cpu_rg=4,cpu_merge=4,cpu_flagstat=2,
                           flow_base_path="/scratch/iacs/iacs_dep/sseth/flows",execute=FALSE){
    ## ----------------- paths and variables
    flowname <- "aln_merge"
    java_tmp <- "java_tmp"
    q_obj <- queue(type=q_type,queue=q_queue)
    ## ----------------- cpu & mem
    flows_desc <- sprintf("%s/%s/aln-merge",project, sample)
    ##junk <- sapply(1:len, function(i){
    if(missing(se.or.pe)) se.or.pe = ifelse(length(fqs2)==length(fqs1), "PE", "SE")
    ## ----------------- setup RG tag variables
    rgid=sprintf("%s_%s", runid, lane, index) ## stays same even if demultplexing is wrong
    rglb=rgsm=sample
    rgpu=lane
    ## ----------------- create paths
    dir.create(bampath, showWarnings=FALSE, recursive=TRUE)
    dir.create(flow_base_path, showWarnings=FALSE, recursive=TRUE)
    ## ----------------- START PIPELINES
    ## ----------------- START alignment
    sai_files1=file.path(gsub(".fastq.gz",".sai",basename(fqs1)))
    sai_files2=file.path( gsub(".fastq.gz",".sai",basename(fqs2)))
    bam_files=file.path(gsub(".fastq.gz",".bam",basename(fqs1)))
    cmd_aln1 <- sprintf("%s aln %s %s %s > %s",bwa_exe, bwa_aln_opts,reflib,fqs1, sai_files1)
    cmd_aln2 <- sprintf("%s aln %s %s %s > %s",bwa_exe, bwa_aln_opts,reflib,fqs2, sai_files2)
    cmd_sampe <- sprintf("%s sampe %s %s %s %s %s %s | %s view -Shu - > %s",
                         bwa_exe,bwa_sampe_opts,reflib,sai_files1,sai_files2,fqs1,fqs2,samtools_exe, bam_files)
    ## ----------------- START read group
    bamrg_files=sprintf("%s_rg.bam",bam_files)
    cmd_fixrg <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/AddOrReplaceReadGroups.jar INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT",
                         java_exe, java_mem, java_tmp, picard_dir, bam_files, bamrg_files, rgid, rglb, platform, rgpu, rgsm,
                         center)
    ## ----------------- START merging
    mergebam <- file.path(out_bampath,out_bam)
    bam_list <- paste("INPUT=", bamrg_files, sep = "", collapse = " ")
    java_mem <- "-Xmx8g";
    cmd_merge <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/MergeSamFiles.jar %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
                         java_exe, java_mem, java_tmp, picard_dir, bam_list, mergebam)
    ## ----------------- samtools flagstat
    flagstat_file <- gsub("bam$","flagstat", mergebam)
    cmd_flagstat <- sprintf("%s flagstat %s > %s",samtools_exe, mergebam, flagstat_file)
    ## ----------------- cleanup
    cmd_cleanup <- "rm -r ../tmp"
    ## ----------------- Define the jobs
    j_obj_aln1 <- job(q_obj=q_obj, name="aln1", cmds=cmd_aln1,  cpu=cpu_aln, submission_type="scatter") ## no dependency
    j_obj_aln2 <- job(q_obj=q_obj, name="aln2", cmds=cmd_aln2,  cpu=cpu_aln, submission_type="scatter") ## no dependency
    j_obj_sampe <- job(q_obj=q_obj, name="sampe", cmds=cmd_sampe,  cpu=cpu_sampe, submission_type="scatter",
                    previous_job=c("aln1", "aln2"), dependency_type="serial")
    j_obj_rg <- job(q_obj=q_obj, name="rg", cmds=cmd_fixrg,  cpu=cpu_rg, submission_type="scatter",
                    previous_job="sampe", dependency_type="serial")
    j_obj_merge <- job(q_obj=q_obj, name="merge", cmds=cmd_merge, cpu=cpu_merge, submission_type="serial",
                       previous_job="rg", dependency_type="gather")
    j_obj_flagstat <- job(q_obj=q_obj, name="flagstat", cmds=cmd_flagstat, cpu=cpu_flagstat,
                          submission_type="serial",
                          previous_job="merge", dependency_type="serial")
    j_obj_cleanup <- job(q_obj=q_obj, name="cleanup", cmds=cmd_cleanup, cpu=cpu_flagstat,
                         submission_type="serial",
                         previous_job="flagstat", dependency_type="serial")
    f_obj <- flow(jobs=list(j_obj_aln1, j_obj_aln2, j_obj_sampe, j_obj_rg, j_obj_merge, j_obj_flagstat, j_obj_cleanup),
                  name=sprintf("%s-%s",flowname, sample),
                  mode="scheduler", flow_base_path=flow_base_path, desc=flows_desc)
    #debug(.submit_flow)
    if(execute) f_obj <- submit_flow(f_obj, execute=TRUE)
    return(f_obj)
}

## passes on arguments to flow
submit_aln_merge <- function(samplematfile, gccpath, out_bampath, se.or.pe,
                             flow_base_path="/scratch/iacs/iacs_dep/sseth/flows",
                             bwa_exe = "/scratch/rists/hpcapps/x86_64/bwa/0.7.9a/bwa",
                             samtools_exe = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools",
                             bwa_aln_opts = "-l 50 -n 6 -k 4 -t 6",
                             bwa_sampe_opts="",
                             java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                             picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                             reflib = "/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/bwa/0.7.5a/Homo_sapiens_assembly19.fasta",
                             q_type="torque",q_queue="iacs",...){
    if(missing(gccpath)) gccpath = "/scratch/iacs/gcc"
    mat <- read.csv(samplematfile, as.is=TRUE)
    samples <- unique(mat$samplename)
    len=length(samples)
    ## ----------------- paths and variables
    q_obj <- queue(type=q_type,queue=q_queue)
    ## ----------------- cpu & mem
    ##junk <- sapply(1:len, function(i){
    tmp <- as.list(match.call(expand.dots=TRUE))[-1]
    tmp <- lapply(tmp,eval, sys.frame(-1)) ## by getting the values from a frame above
    params <- formals(flow_aln_merge)
    tmp <- tmp[names(tmp) %in% names(params)]
    params[names(tmp)]=tmp             #replace in formals
    ret <- list()
    for(i in 1:len){
        cat(".")
        ## ------------------  DEFINE PATHS
        sample <- samples[i]
        sample.mat <- subset(mat,samplename==sample)
        s.runid <- sample.mat$runid[1];s.lane=sample.mat$lane[1]
        s.project=sample.mat$project[1];s.subproject=sample.mat$subproject[1]
        fqs1=as.c(subset(sample.mat,read==1)$files)
        fqs2=as.c(subset(sample.mat,read==2)$files)
        sorted_bam <- sample.mat$sorted_bam[1]
        params$fqs1=fqs1;params$fqs2=fqs2;params$out_bam=sorted_bam
        params$sample=sample;params$project=s.project
        f_obj <- do.call(flow_aln_merge, args = params)
        ret <- c(ret, f_obj) ## list of flows
    }
    return(ret)
}


if(FALSE){

    source("~/Dropbox/public/github.flow/R/class-def.R")
    source("~/Dropbox/public/github.flow/R/class-funcs.R")
    source("~/Dropbox/iacsSVN/RPacks/filenames/R/fastq_files.R")
    require(uuid)
    ## debug(.submit_flow)
    ## debug(.submit_job)
    ## debug(.create_queue_cmd)
    f_obj <- flow_aln_merge(fqs1=fqs1, fqs2=fqs2, out_bam=out_bam, sample="1-NS", project="MTSCC",
                            q_type="lsf",q_queue="normal",execute=TRUE,
                            out_bampath="/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_cSCC_Ni56-2/bams",
                            bwa_aln_opts = "-l 50 -n 6 -k 4 -t 16",
                            cpu_aln=16,cpu_sampe=1,cpu_rg=1,cpu_merge=1,cpu_flagstat=1)
    ## out <- submit_flow(f_obj)
    ## if(execute) f_obj <- submit_flow(f_obj, execute=TRUE)


    se.or.pe <- "PE"
    species <- "human"
    runid <- "mt_scc"
    ## undebug(create_sample_mat)
    ## mat <- create_sample_mat(path="/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-NI56-1/orig_fqs", project="MTSCC",
    ##                   subproject="NI56A", outpath="~/projects/analysis_ss.mt.scc/files")
    ## mat <- filenames::create_sample_mat(path="/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_cSCC_Ni56-2/orig_fqs",
    ##                                     project="MTSCC",subproject="NI562", outpath="~/projects/analysis_ss.mt.scc/files")
    ## mat <- filenames::create_sample_mat(path="/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-prog_HC3/orig_fqs",
    ##                                     project="MTSCC",subproject="HC3", outpath="~/projects/analysis_ss.mt.scc/files")

    samplematfile <- "~/projects/analysis_ss.mt.scc/files/MTSCC_NI56A_Project_VC_SCC-NI56-1_sample_mat.csv"
    bampath <- "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-NI56-1/bams"

    samplematfile <- "~/projects/analysis_ss.mt.scc/files/MTSCC_NI562_Project_VC_cSCC_Ni56-2_sample_mat.csv"
    bampath <- "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_cSCC_Ni56-2/bams"

    ## samplematfile <- "~/projects/analysis_ss.mt.scc/files/MTSCC_HC3_Project_VC_SCC-prog_HC3_sample_mat.csv"
    ## bampath <- "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-prog_HC3/bams"

    source("~/Dropbox/public/github.r-ngs-utils/R/pipe_fastq_bwa.R")
    ## debug(submit_aln_merge)
    f_objs <- submit_aln_merge(samplematfile, out_bampath=bampath,
                     q_type="lsf",q_queue="normal",execute=TRUE,
                     bwa_aln_opts = "-l 50 -n 6 -k 4 -t 4",
                     cpu_aln=4,cpu_sampe=1,cpu_rg=1,cpu_merge=1,cpu_flagstat=1)

}



##
if(FALSE){

    samplematfile <- "~/projects/analysis_ss.mt.scc/files/MTSCC_NI562_Project_VC_cSCC_Ni56-2_sample_mat.csv"
    bampath <- "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_cSCC_Ni56-2/bams"
    i=1
    mat <- read.csv(samplematfile, as.is=TRUE)
    samples <- unique(mat$samplename)
    sample <- samples[i]
    sample.mat <- subset(mat,samplename==sample)
    runid <- sample.mat$runid[1];lane=sample.mat$lane[1]
    project=sample.mat$project[1];subproject=sample.mat$subproject[1]
    fqs1=as.c(subset(sample.mat,read==1)$files)
    fqs2=as.c(subset(sample.mat,read==2)$files)
    out_bam <- sample.mat$sorted_bam[1]


    bwa_exe <- "/scratch/rists/hpcapps/x86_64/bwa-0.7.5a"
    samtools_exe <- "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools"
    bwa_aln_opts="-l 50 -n 6 -k 4 -t 6"
    bwa_sampe_opts=""
    java_exe <- "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java"
    tmppath <- "tmp"
    picard_dir <- "/scratch/rists/hpcapps/x86_64/picard/1.112"
    reflib <- "/scratch/rists/hpcapps/reference/human/b37/indexes/bwa/0.7.5a/human_g1k_mt_NC_012920.1.fasta"
    rgpl=platform='illumina'
    rgcn=center='IACS'
    flow_base_path="/scratch/iacs/iacs_dep/sseth/flows"
    q_type="lsf"; q_queue="normal"

}


