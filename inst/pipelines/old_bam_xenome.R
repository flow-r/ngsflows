
if(FALSE){

    require(flow)
    q_obj <- queue(type="lsf",queue="normal")
    bam_file="/IACS1/GCC/LevelII/130820_SN1440_0162_BC29KLACXX/GIULIO-GSC-ES-811-10-01D_130820_SN1440_0162_BC29KLACXX_s_8_TCTTCA.rg.sorted.recalibed.bam"
    sample="GIULIO-GSC-ES-811-10-01D";project="project"
    java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java"
    java_mem_bam_fq= "-Xmx4g"
    java_mem_xenome = "-Xmx12g"
    java_tmp = "java_tmp"
    picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112"
    xenome_exe = "/scratch/rists/hpcapps/x86_64/xenome/1.0.1-r/xenome"
    xen_index_broad_hg19_mm9 = "/IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index"
    xen_fastq_exe = "/IACS1/home/sseth/projects/tools_xenome/xenome.to.fastq.sh"
    xen_filter_bam_exe = "/IACS1/home/sseth/Dropbox/projects/tools_xenome/filter.xeno.mmuReads.sseth.py"
    cpu_bam_fq = 1; cpu_xenome = 24; cpu_xen_fq = 1; cpu_filter = 1
    flow_base_path="/scratch/iacs/iacs_dep/sseth/flows"

    fobj <- flow_bam_xenome(bam_file=bam_file, q_obj=q_obj)
    plot(fobj, pdf = TRUE, detailed = TRUE)


    require(flow)
    source("~/projects/tools_flowpipes/pipe_bam_xenome.R")
    source("~/projects/tools_flowpipes/samplesheet_funcs.R")
    source("~/Dropbox/public/github_flow/R/class-funcs.R")
    source("~/Dropbox/public/github_flow/R/class-def.R")
    #sheet <- "/scratch/iacs/iacs_dep/sseth/flows/GIULIO-GSC-JG/paired_sample_sheet.tsv"
    sheet <- "/scratch/iacs/iacs_dep/sseth/flows/LYNDA-FORD/paired_sample_sheet.tsv"
    flows <- submit_bam_xenome(sheet, execute=TRUE, q_type="torque",q_queue="long")

}
#########------------ HPCC; long/iacs
## cpu_bam_fq = 4, cpu_xenome = 24, cpu_xen_fq = 1, cpu_filter = 1,
#########------------ LSF
## cpu_bam_fq = 1, cpu_xenome = 16, cpu_xen_fq = 1, cpu_filter = 1,
flow_bam_xenome <- function(bam_file, q_obj=q_obj,
                            project="project",sample="sample",
                            flow_base_path="/scratch/iacs/iacs_dep/sseth/flows",
                            java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                            java_mem_bam_fq= "-Xmx4g",
                            java_mem_xenome = "-Xmx12g",
                            java_tmp = "java_tmp",
                            picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                            xenome_exe = "/scratch/rists/hpcapps/x86_64/xenome/1.0.1-r/xenome",
                            xen_index_broad_hg19_mm9 = "/IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index",
                            xen_fastq_exe = "~/projects/tools_xenome/xenome.to.fastq.sh",
                            xen_filter_bam_exe = "~/projects/tools_xenome/filter.xeno.mmuReads.sseth.py",
                            cpu_bam_fq = 1, cpu_xenome = 24, cpu_xen_fq = 1, cpu_filter = 1,
                            execute=FALSE){
    ## ------- inputs
    flow_name = "bam_xenome"
    fq1="read1.fq";fq2="read2.fq"
    ## ------- outputs
    out_bam_file <- gsub("bam$", "xeno_bam", bam_file)
    ## ------- convert bam to fastq
    tool="SamToFastq.jar"
    java_params="-XX:-UseGCOverheadLimit -XX:-UseParallelGC"
    cmd_bam_fq <- sprintf("time %s %s %s -Djava.io.tmpdir=%s -jar %s/%s INPUT=%s FASTQ=%s SECOND_END_FASTQ=%s VALIDATION_STRINGENCY=LENIENT",
                          java_exe, java_mem_bam_fq, java_params, java_tmp, picard_dir, tool, bam_file, fq1, fq2)
    ## ------- run xenome on the fastq files
    xen_prefix = sample;#gsub(".bam$","",basename(bam_file))
    cmd_xenome <- sprintf("%s classify --num-threads %s -P %s --pairs --graft-name hg19 --host-name mm9 --output-filename-prefix %s -i %s -i %s -v ",
                          xenome_exe, cpu_xenome, xen_index_broad_hg19_mm9, xen_prefix, fq1, fq2)
    ## ------- convert xen fastq to be valid
    xen.out.suffix <- c("mm9_1.fastq","mm9_2.fastq","neither_1.fastq","neither_2.fastq",
                        "ambiguous_1.fastq","ambiguous_2.fastq","both_1.fastq","both_2.fastq")
    cmd_xen_fq <- sprintf("%s %s > %s", xen_fastq_exe, xen.out.suffix,
                          gsub("(.*)_([1-2]{1}).fastq","\\1_end\\2.fq", xen.out.suffix))
    ## ------- filter the bam file using the xen fastq
    env="export PYTHONPATH=$PYTHONPATH:/IACS1/lib/python2.7/site-packages/"
    ## -----------------------------------------env   python script     outpath  prefix    bamfile    output_bam
    cmd_filter <- sprintf("%s; %s %s %s %s %s", env, xen_filter_bam_exe, "./", xen_prefix, bam_file, out_bam_file)
    ## ------- job
    j_obj_bam_fq <- job(q_obj=q_obj, cmds = cmd_bam_fq, cpu = cpu_bam_fq, submission_type="serial", name="bam_fq")
    j_obj_xenome <- job(q_obj = q_obj, cmds = cmd_xenome, cpu = cpu_xenome, submission_type = "serial",
                        name = "xenome", dependency_type = "serial", previous_job = "bam_fq")
    j_obj_xen_fq <- job(q_obj = q_obj, cmds = cmd_xen_fq, cpu = cpu_xen_fq, submission_type = "scatter",
                        name = "xen_fq", dependency_type = "burst", previous_job = "xenome")
    j_obj_filter <- job(q_obj = q_obj, cmds = cmd_filter, cpu = cpu_filter, submission_type = "serial",
                        name = "filter", dependency_type = "gather", previous_job = "xen_fq")
    ## ------- flow
    flow_desc <- sprintf("%s/%s/%s",project, sample, flow_name)
    f_obj <- flow(jobs=list(j_obj_bam_fq, j_obj_xenome, j_obj_xen_fq, j_obj_filter), name = flow_name,
                  desc=flow_desc, mode="scheduler", flow_base_path=flow_base_path)
    f_obj <- submit_flow(f_obj, execute=execute)
    return(f_obj)
}

submit_bam_xenome <- function(paired_samplesheet,
                              flow_base_path="/scratch/iacs/iacs_dep/sseth/flows",
                              java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                              java_mem_bam_fq= "-Xmx4g",
                              java_mem_xenome = "-Xmx12g",
                              picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                              xenome_exe = "/scratch/rists/hpcapps/x86_64/xenome/1.0.1-r/xenome",
                              xen_index_broad_hg19_mm9 = "/IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index",
                              xen_fastq_exe = "~/projects/tools_xenome/xenome.to.fastq.sh",
                              xen_filter_bam_exe = "~/projects/tools_xenome/filter.xeno.mmuReads.sseth.py",
                              cpu_bam_fq = 1, cpu_xenome = 24, cmd_xen_fq = 1, cmd_filter = 1,
                              q_type="torque",q_queue="iacs",...){
    mat <- read_paired_samplesheet(paired_samplesheet)$single_mat
    samples <- unique(c(mat$samplename, mat$refname))
    len=length(samples)
    ## ----------------- paths and variables
    q_obj <- queue(type=q_type,queue=q_queue)
    ## ----------------- cpu & mem
    ##junk <- sapply(1:len, function(i){
    tmp <- as.list(match.call(expand.dots=TRUE))[-1]
    tmp <- lapply(tmp,eval, sys.frame(-1)) ## by getting the values from a frame above
    params <- formals(flow_bam_xenome)
    tmp <- tmp[names(tmp) %in% names(params)]
    params[names(tmp)]=tmp             #replace in formals
    ret <- list()
    for(i in 1:len){
        cat(".")
        ## ------------------  DEFINE COMMON PATHS
        sample <- samples[i]
        sample.mat <- subset(mat,samplename==sample)
        s.project <- sample.mat$project[1]
        s.bam <- sample.mat$bam[1]
        ## ------------------  DEFINE COMMON PATHS into PARAMS
        params$bam_file=s.bam;params$sample=sample;params$project=s.project
        params$q_obj=q_obj
        f_obj <- do.call(flow_bam_xenome, args = params)
        ret <- c(ret, f_obj) ## list of flows
    }
    return(ret)
}
