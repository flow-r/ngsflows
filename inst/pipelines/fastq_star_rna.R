## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com

## =============== RNASeq pipeline which uses STAR and then unified genotyper for SNP calling
## --------------- setup for MDA HPCC
## cpu_fixrg = 4

rna_star_preprocess <- function(fastq1,fastq2,out_prefix="sample",
                              project = "project",sample = "sample",
                              star_exe = "/scratch/rists/hpcapps/x86_64/STAR/2.3.0e/STAR",
                              reference_star = "/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/STAR/2.3.0e",
                              reffa = "/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta",
                              q_type="lsf", q_queue="normal",
                              picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                              java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                              java_mem_rg = "-Xmx2g",
                              cpu_merge=2, cpu_star=16, cpu_fixrg = 4,
                              ## RG tags stuff
                              platform  = 'illumina',center = 'IACS',runid = "runid",lane = "0", index = "index"){
    ## ----------------- paths and variables
    flowname <- "rna_star"
    java_tmp <- "java_tmp"
    q_obj <- queue(type=q_type,queue=q_queue)
    flow_desc <- sprintf("%s/%s/%s",project, sample, flowname)
    ## ----------------- setup RG tag variables
    rgid=sprintf("%s_%s", runid, lane, index)
    ## stays same even if demultplexing is wrong
    rglb=rgsm=sample;rgpu=lane
    ## ----------------- merge fastq files
    #cmd_merge <- merge_fastqs(run=runid, proj=project, sample=sample)
    ## ----------------- STAR
    outdir = "."
    tmp <- star_pipe(fastq1=fastq1, fastq2=fastq2, star_exe=star_exe,
                     reffa=reffa, genome_dir=reference_star, threads=cpu_star)
    cmd_star <- tmp$cmds
    ## ----------------- fix RG commands
    sam=tmp$outfile;
    sorted_bam=sprintf("%s.sorted.bam",out_prefix)
    cmd_fixrg <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/AddOrReplaceReadGroups.jar INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT",
                         java_exe, java_mem_rg, java_tmp, picard_dir, sam, out_prefix, rgid, rglb, platform, rgpu, rgsm,
                         center)
    ## ----------------- preprocess
    recal_bam=sprintf("%s.sorted.recalibed.bam",out_prefix)
    tmp <- bam_preprocess(inbam=sorted_bam, outbam=recal_bam, q_obj=q_obj)
    jobs_preprocess <- tmp$flow@jobs
    jobs_preprocess$markdup@previous_job="fixrg"
    jobs_preprocess$markdup@dependency_type="serial"
    ## ---------------- make flow
    ## j_obj_merge <- job(q_obj=q_obj, name="merge", cmds=cmd_merge,  cpu=cpu_merge,
    ##                    submission_type="serial")
    j_obj_star <- job(q_obj=q_obj, name="star", cmds=cmd_star,  cpu=cpu_star,
                      submission_type="serial")
    j_obj_fixrg <- job(q_obj=q_obj, name="fixrg", cmds=cmd_fixrg,  cpu=cpu_fixrg,
                       previous_job="star", dependency_type="serial")
    jobs=c(j_obj_star, j_obj_fixrg, unlist(jobs_preprocess))
    fobj <- flow(jobs=jobs,name=flowname, mode="scheduler", desc=flow_desc)
    return(fobj)
}

if(FALSE){

    fastq1='/scratch/iacs/gcc/leveli/140517_SN1222_0251_AC3LTPACXX/Project_IACS-Tim/Sample_TH-AC10-60-01R/TH-AC10-60-01R_GTGGCC_L008_R1_001.fastq.gz'
    fastq2='/scratch/iacs/gcc/leveli/140517_SN1222_0251_AC3LTPACXX/Project_IACS-Tim/Sample_TH-AC10-60-01R/TH-AC10-60-01R_GTGGCC_L008_R2_001.fastq.gz'

    require(flow)
    source("~/iacsSVN/RScripts/Rfuncs/star_func.R")
    source("~/Dropbox/projects/tools_flowpipes/pipe_rna_star_variants.R");
    source("~/Dropbox/public/github_r-ngs-utils/R/funcs_bam_preprocess.R")
    ## debug(rna_star_variants)
    fastq1="/scratch/iacs/gcc/levelimergedfqs/140530_SN1222_0256_AC4DMYACXX/IACS-Tim-TH-AC01-60-01R_140530_SN1222_0256_AC4DMYACXX_s_5_GCCAAT.R1.fastq"
    fastq2="/scratch/iacs/gcc/levelimergedfqs/140530_SN1222_0256_AC4DMYACXX/IACS-Tim-TH-AC01-60-01R_140530_SN1222_0256_AC4DMYACXX_s_5_GCCAAT.R2.fastq"
    fobj=rna_star_preprocess(fastq1=fastq1,fastq2=fastq2, project='ALE',
        sample='Sample_TH-AC10-60-01R', runid="140530_SN1222_0256_AC4DMYACXX")
    fobj <- submit_flow(fobj, execute=TRUE)


    star_exe = "/scratch/rists/hpcapps/x86_64/STAR/2.3.0e/STAR"
    reference_star = "/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/STAR/2.3.0e"
    cpu_star = 16
    reffa = "/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta"
    fastq1 = "fastq1"
    fastq2 = "fastq2"
    outdir = "."

    project = "project"
    sample = "sample"
    q_type="torque"
    q_queue="iacs"

    sorted_bam = "sample.sorted.bam"

    runid = "runid"
    lane = "0"
    index = "index"
    sample = "sample"
    platform  = 'illumina'
    center = 'IACS'

    require(flow)
    source("~/Dropbox/public/github.r-ngs-utils/R/funcs_bam_preprocess.R")
    bam="/scratch/iacs/gcc/levelii/140603_SN208_0516_AC4HL2ACXX/ANDY-Futreal-AF-780702-01-01D_140603_SN208_0516_AC4HL2ACXX_s_4_AACGCTTA.rg.sorted.bam"
    outbam=gsub(".bam$",".reacalibed.bam", bam)
    q_obj <- queue(type="torque",queue="iacs")

    debug(bam_preprocess)
    tmp <- bam_preprocess(inbam=bam, outbam=outbam, q_obj=q_obj)

}
