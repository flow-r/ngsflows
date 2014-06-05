## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------


pipe_aln_merge <- function(samplematfile, gccpath, bampath, qapath, gatkpath, se.or.pe,
                           flow_base_path="/scratch/iacs/iacs_dep/sseth/flows",
                           bwa_exe = "/scratch/iacs/tmp/bwa-0.7.9a/bwa",
                           samtools_exe = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools",
                           bwa_aln_opts = "-l 50 -n 6 -k 4 -t 6",
                           bwa_sampe_opts="",
                           java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                           picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                           reflib = "/scratch/rists/hpcapps/reference/human/b37/indexes/bwa/0.7.5a/human_g1k_mt_NC_012920.1.fasta",
                           platform  = 'illumina',
                           center = 'IACS',
                           tmppath = "tmp",
                           q_type="torque",q_queue="iacs",
                           execute=FALSE){
    if(missing(gccpath)) gccpath = "/scratch/iacs/gcc"
    mat <- read.csv(samplematfile, as.is=TRUE)
    samples <- unique(mat$samplename)
    len=length(samples)
    ## ----------------- paths and variables
    flowname <- "aln_merge"
    q_obj <- queue(type=q_type,queue=q_queue)
    ## ----------------- cpu & mem
    java_mem <- "-Xmx4g"
    cpu_aln <- 12
    cpu_rg <- 4
    cpu_merge <- 12
    cpu_flagstat <- 2
    ##junk <- sapply(1:len, function(i){
    ret <- list()
    for(i in 1:len){
        ## ------------------  DEFINE PATHS
        sample <- samples[i]
        sample.mat <- subset(mat,samplename==sample)
        s.runid <- sample.mat$runid[1];s.lane=sample.mat$lane[1]
        s.project=sample.mat$project[1];s.subproject=sample.mat$subproject[1]
        if(missing(bampath)) bampath = file.path(gccpath,"levelii",s.runid) ## folder for samplemat
        if(missing(qapath)) qapath = file.path(gccpath,"qa",s.runid,sample)
        if(missing(gatkpath)) gatkpath = file.path(gccpath,"leveliii",s.runid, "gatk.UniTyper") ## folder for samplemat
        if(missing(se.or.pe)) se.or.pe = ifelse(length(unique(sample.mat$read))==2, "PE", "SE")
        out_basename <- as.c(sample.mat$out_basename[1])
        sorted_bam <- sample.mat$sorted_bam[1]
        fqs1=as.c(subset(sample.mat,read==1)$files)
        fqs2=as.c(subset(sample.mat,read==2)$files)
        ## ----------------- setup variables
        cat("Sample:", sample,"\n")
        s.flow_base_path <- sprintf("%s/%s/%s",flow_base_path,s.project,sample)
        rgid=rglb=rgsm=sample
        rgpu=s.lane
        ## ----------------- create paths
        dir.create(bampath, showWarnings=FALSE, recursive=TRUE)
        dir.create(qapath, showWarnings=FALSE, recursive=TRUE)
        dir.create(gatkpath, showWarnings=FALSE, recursive=TRUE)
        dir.create(s.flow_base_path, showWarnings=FALSE, recursive=TRUE)
        ## ----------------- START PIPELINES
        ## ----------------- START alignment
        out_files=file.path(tmppath, gsub(".fastq.gz",".bam",basename(fqs1)))
        cmd_aln1 <- sprintf("%s aln %s %s %s",bwa_exe, bwa_aln_opts,reflib,fqs1)
        cmd_aln2 <- sprintf("%s aln %s %s %s",bwa_exe, bwa_aln_opts,reflib,fqs2)
        cmd_aln <- sprintf("%s sampe %s %s <(%s) <(%s) %s %s | %s view -Shu - > %s",
                           bwa_exe,bwa_sampe_opts,reflib,cmd_aln1,cmd_aln2,fqs1,fqs2,samtools_exe, out_files)
        in_files=out_files ### for the next step
        ## ----------------- START read group
        out_files=sprintf("%s_rg.bam",out_files)
        cmd_fixrg <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/AddOrReplaceReadGroups.jar INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT",
                             java_exe, java_mem, tmppath, picard_dir, in_files, out_files, rgid, rglb, platform, rgpu, rgsm,
                             center)
        ## ----------------- START merging
        in_files=out_files; out_file <- file.path(bampath,sorted_bam)
        bam_list <- paste("INPUT=", in_files, sep = "", collapse = " ")
        java_mem <- "-Xmx8g";
        cmd_merge <- sprintf("%s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/MergeSamFiles.jar %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
                             java_exe, java_mem, picard_dir, bam_list, out_file)
        ## ----------------- samtools flagstat
        in_file=out_file;out_file <- gsub("bam$","flagstat", in_file)
        cmd_flagstat <- sprintf("%s flagstat %s > %s",samtools_exe, in_file, out_file)
        ## ----------------- Define the jobs
        j_obj_aln <- job(q_obj=q_obj,cmds=cmd_aln, cpu=cpu_aln, submission_type="scatter", name="aln") ## no dependency
        j_obj_rg <- job(q_obj=q_obj,cmds=cmd_fixrg, cpu=cpu_rg, submission_type="scatter",dependency_type="serial", name="rg",
                        previous_job="aln")
        j_obj_merge <- job(q_obj=q_obj, cmds=cmd_merge, submission_type="serial",dependency_type="gather",
                           name="merge", cpu=cpu_merge, previous_job="rg")
        j_obj_flagstat <- job(q_obj=q_obj,cmds=cmd_flagstat, submission_type="serial",dependency_type="serial",
                           name="flagstat", cpu=cpu_flagstat, previous_job="merge")
        f_obj <- flow(jobs=list(j_obj_aln, j_obj_rg, j_obj_merge, j_obj_flagstat),
                      name=sprintf("%s-%s",flowname, sample),
                      mode="scheduler", flow_base_path=s.flow_base_path)
        if(execute) f_obj <- submit_flow(f_obj, execute=TRUE)
        ret <- c(ret, f_obj) ## list of flows
    }
    return(ret)
}

as.c=as.character

if(FALSE){

    source("~/Dropbox/public/github.flow/R/class-def.R")
    source("~/Dropbox/public/github.flow/R/class-funcs.R")
    source("~/Dropbox/iacsSVN/RPacks/filenames/R/fastq_files.R")
    source("~/Dropbox/public/github.r-ngs-utils/R/pipe_fastq_bwa.R")
    require(uuid)

    ## bwa_exe <- "/scratch/rists/hpcapps/x86_64/bwa-0.7.5a"
    ## samtools_exe <- "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools"
    ## bwa_aln_opts="-l 50 -n 6 -k 4 -t 6"
    ## bwa_sampe_opts=""
    ## java_exe <- "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java"
    ## tmppath <- "tmp"
    ## picard_dir <- "/scratch/rists/hpcapps/x86_64/picard/1.112"
    ## reflib <- "/scratch/rists/hpcapps/reference/human/b37/indexes/bwa/0.7.5a/human_g1k_mt_NC_012920.1.fasta"
    ## rgpl=platform='illumina'
    ## rgcn=center='IACS'

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

    samplematfile <- "~/projects/analysis_ss.mt.scc/files/MTSCC_HC3_Project_VC_SCC-prog_HC3_sample_mat.csv"
    bampath <- "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-prog_HC3/bams"

    source("~/Dropbox/public/github.r-ngs-utils/R/pipe_fastq_bwa.R")
    source("~/Dropbox/public/github.flow/R/class-def.R")
    source("~/Dropbox/public/github.flow/R/class-funcs.R")
    source("~/Dropbox/public/github.flow/R/generic-funcs.R")
    require(uuid)
    ##f_obj <- pipe_aln_merge(samplematfile, bampath=bampath, execute=FALSE)
    f_objs <- pipe_aln_merge(samplematfile, bampath=bampath, execute=TRUE, q_type="lsf",q_queue="normal")

}


