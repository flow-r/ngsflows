## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------

## barcodes: CGAA
## merge.dat splot by samples
get_trim_aln_cmds <- function(samplematfile,trim.table,pe,reflibs,
                              samtools_exe="/IACS1/NGS/samtools/samtools",
                              barcodes='/IACS1/home/sseth/projects/analysis_ac.seq.tag/input.files/barcodes.txt',
                              splitter='/IACS1/home/sseth/code/perl/fastx_barcode_splitter.pl',
                              bowtie="/IACS1/apps/bowtie2-2.1.0/bowtie2",
                              picard_dir="/scratch/rists/hpcapps/x86_64/sequence/picard-tools-1.48",
                              aln_options="-N 1 -L 18 --end-to-end",
                              java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                              execute=FALSE, ## whether or not execute the flow
                              tmppath="tmp",## relative path to CWD which temp files are kept
                              filter.cores=24,trim.cores=24,bowtie.cores=24){
    mat <- read.csv(samplematfile, as.is=TRUE)
    samples <- unique(mat$samplename)
    len <- length(samples)
    ## q_obj <- queue(type="torque",queue="iacs")
    flowname <- "shrna"
    cpu_filtered <- 4
    cpu_aln <- 12
    cpu_uniq <- 4
    cpu_merge <- 12
    cpu_flag <- 2
    cpu_trim <- 4
    cpu_idx <- 4
    ret <- list()
    for(i in 1:len){
        ## junk <- lapply(1:length(merge.dat),function(i){
        sample <- samples[i]
        sample.mat <- subset(mat,samplename==sample)
        s.runid <- sample.mat$runid[1];s.lane=sample.mat$lane[1]
        s.project=sample.mat$project[1];s.subproject=sample.mat$subproject[1]
        flow_base_path <- sprintf("~/tmp/flows/%s/%s",s.project,sample)
        fqs1=as.c(subset(sample.mat,read==1)$files)
        filtered_files=file.path(tmppath, gsub(".fastq.gz",".pass.fastq",
            basename(fqs1)))
        trimmed_files1=file.path(tmppath, gsub(".fastq.gz",".trimmed1.fastq",
            basename(fqs1)))
        trimmed_files2=file.path(tmppath, gsub(".fastq.gz",".trimmed2.fastq",
            basename(fqs1)))
        cmd_filter=sprintf("gunzip -c %s | %s --bcfile %s --prefix %s --suffix .fastq --barcode_start 18 --mismatches 1",
            fqs1,splitter,barcodes,gsub("pass.fastq","",filtered_files))
        ## trim them:
        cmd_trim1 <- unlist(trim.fastq(input.files=filtered_files,
                                       output.files=trimmed_files1,
                                       from=trim.table[1,1],to=trim.table[1,2],
                                       qual="Q33",cores=24,execute=FALSE))
        cmd_trim2 <- unlist(trim.fastq(input.files=filtered_files,
                                       output.files=trimmed_files2,
                                       from=trim.table[3,1],to=trim.table[3,2],
                                       qual="Q33",cores=24,execute=FALSE))
        ## Align them:
        bam_files=file.path(tmppath, gsub(".fastq.gz",".bam",basename(fqs1)))
        reflib=reflibs[i]
        if(pe[i]){
            cmd_aln <- sprintf("%s --time %s -x %s -1 %s -2 %s -p %s | %s view -Shu - | %s sort - %s",
                               bowtie,aln_options,reflib,trimmed_files1,
                               trimmed_files2,bowtie.cores,samtools_exe, samtools_exe,
                               gsub(".bam$","",bam_files))
        }else{
            cmd_aln <- sprintf("%s --time %s -x %s -U %s -p %s | %s view -Shu - | %s sort - %s",
                               bowtie,aln_options,reflib,trimmed_files2,
                               bowtie.cores,samtools_exe,samtools_exe,
                               gsub(".bam$","",bam_files))
        }
        ## ------- Merge bams
        ## ------- inputs
        bam_list <- paste("INPUT=", bam_files, sep = "", collapse = " ")
        ## ------- outputs
        sorted_bam <- sample.mat$sorted_bam[1]
        cmd_merge <- sprintf("%s %s -Djava.io.tmpdir=tmp -jar %s/MergeSamFiles.jar %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
                             java_exe, java_mem, picard_dir, bam_list, sorted_bam)
        ## Get counts
        out_counts <- gsub(".bam$",".counts.txt",sorted_bam)
        ## use only uniquely aligned
        cmd_uniq <- sprintf("%s view -f 0x0040 %s | awk '{clones[$3\"\t\"$7]++}END{for (clone in clones) print clone,\"\t\",clones[clone]}' > %s",
                            samtools_exe, sorted_bam,out_counts)
        idxstats <- gsub(".bam$",".idxstats.txt",sorted_bam)
        cmd_idx <- sprintf("samtools idxstats %s > %s",sorted_bam,idxstats)
        flagstat <- gsub(".bam$",".flagstat.txt",sorted_bam)
        cmd_flag <- sprintf("samtools flagstat %s > %s",sorted_bam,flagstat)
        ## ----------------- Define the jobs
        cmd_trim <- paste(cmd_trim1, cmd_trim2, sep=";")
        q_obj <- queue(type="lsf",queue="normal")
        j_obj_filter <- job(q_obj=q_obj,cmds=cmd_filter,cpu=cpu_filtered,
                            submission_type="scatter", name="filter")
        j_obj_trim <- job(q_obj=q_obj,cmds=cmd_trim,cpu=cpu_trim,
                           name="trim", submission_type="scatter", dependency_type="serial",
                           previous_job="filter")
        j_obj_aln <- job(q_obj=q_obj,cmds=cmd_aln,cpu=cpu_aln, dependency_type="serial",submission_type="scatter", name="aln",
                         previous_job="trim") ## no dependency
        j_obj_merge <- job(q_obj=q_obj, cmds=cmd_merge, submission_type="serial",
                           dependency_type="gather",name="merge", cpu=cpu_merge,
                           previous_job="aln")
        j_obj_uniq <- job(q_obj=q_obj,cmds=cmd_uniq, submission_type="serial",
                          dependency_type="serial",name="uniq", cpu=cpu_uniq,
                          previous_job="merge")
        j_obj_flag <- job(q_obj=q_obj,cmds=cmd_flag,
                          submission_type="serial",dependency_type="serial",
                          name="flagstat", cpu=cpu_flag,
                          previous_job="merge")
        j_obj_idx <- job(q_obj=q_obj,cmds=cmd_idx,
                         submission_type="serial",dependency_type="serial",
                         name="idx", cpu=cpu_idx,
                         previous_job="merge")
        f_obj <- flow(jobs=list(j_obj_filter, j_obj_trim, j_obj_aln,
                          j_obj_merge, j_obj_uniq, j_obj_flag, j_obj_idx),
                      name=sprintf("%s-%s",flowname, sample),
                      mode="scheduler", flow_base_path=flow_base_path)
        names(f_obj@jobs) <- sapply(f_obj@jobs, slot, "name")
        if(execute) f_obj_uuid <- .submit_flow(f_obj, execute=TRUE)
        #f_obj <- .submit_flow(f_obj, execute=TRUE, attach_uuid=FALSE)
        ret <- c(ret, f_obj) ## list of flows
    } ## for loop
    return(ret)
}


if(FALSE){

    source("~/Dropbox/public/github.flow/R/generic-funcs.R")
    source("~/Dropbox/public/github.flow/R/class-def.R")
    source("~/Dropbox/public/github.flow/R/class-funcs.R")
    source("~/Dropbox/iacsSVN/RPacks/filenames/R/fastq_files.R")
    require(mypack);require(uuid)
    base_path <- "/scratch/iacs/iacs_dep/sseth/data/acarugo.read.counting"
    runid <- "Project_AC_evol-GBM_PM28"
    format <- "$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"
    samplematfile <- file.path(base_path, runid,"shrna_13k_Project_AC_evol-GBM_PM28_sample_mat.csv")
    ## mat <- create_sample_mat(path=file.path(base_path,runid,"orig_fqs"),
    ##                          project="shrna",subproject="13k",
    ##                          format=format,
    ##                          outpath=file.path(base_path,runid))
    trim.table=matrix(c(1,18,19,22,23,40),nrow=3,byrow=TRUE)
    scramble <- "/scratch/iacs/iacs_dep/sseth/data/acarugo.read.counting/nontargeting-barcodes-13K-library/nontargeting-barcodes-13K-library-bwt2"
    reflibs=c(rep(scramble,8));pe=rep(TRUE,8)
    ## barcodes='/scratch/iacs/iacs_dep/sseth/projects/analysis_ac.seq.tag/input.files/barcodes.txt'
    ## bowtie="/scratch/rists/hpcapps/x86_64/bowtie2/2.2.2/bowtie2"
    ## picard_dir="/scratch/rists/hpcapps/x86_64/sequence/picard-tools-1.48"
    ## aln_options=" --ff  -N 1 -L 18 --end-to-end"
    ## bowtie.cores <- 16
    ## splitter='/scratch/iacs/iacs_dep/sseth/code/perl/fastx_barcode_splitter.pl'
    ## samtools_exe="/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools"
    source("~/Dropbox/public/github.r-ngs-utils/R/pipe_trim_aln_merge.R")
    f_objs <- get_trim_aln_cmds(samplematfile=samplematfile,trim.table=trim.table, reflibs=reflibs,pe=pe, tmppath="tmp")
    f_obj <- f_objs[[1]]
    plot(f_obj, detailed=TRUE)
}
