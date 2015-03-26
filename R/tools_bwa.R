## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------

require(flow)
setClass("bwa", contains = "job",
         representation(fastq1 = "character", ## submit job
                        fastq2 = "character", ## type of queue
                        paired_end = "logical", 
                        bwa_exe = "character",
                        bwa_command = "character",
                        bwa_ref = "character",
                        bwa_opt = "character"
         )) ## address of head node


## this should accept all the commands and create a command string for job_cmd
bwa <- function(fastq1 = '', fastq2 = '', bwa_exe = 'bwa', bwa_command = c("mem", "aln_sam"), paired_end = TRUE, bwa_opt='',
                bwa_ref = '',...){
  ## other arguments passed on to job class
  bwa_command <- match.arg(bwa_command)
  if(bwa_command == "mem" & paired_end){
    cmds <- sprintf("%s mem %s %s %s %s",
                    bwa_exe, bwa_opt, bwa_ref, fastq1, fastq2)
  }else if(bwa_command == "mem" & !paired_end){
    cmds <- sprintf("%s mem %s %s",
                    bwa_exe, bwa_opt, bwa_ref, fastq1)
  }
  object <- new("bwa", bwa_command=bwa_command, cmds=cmds, name = "bwa", ...)
  return(object)
}

if(FALSE){
  bwa()@cmds
  #cmd_aln1 <- paste(bwapath, "/bwa aln ",bwa_aln_opts," ",reflib," ",fqs1,sep = "")
  
}