## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com

## =============== RNASeq pipeline which uses STAR and then unified genotyper for SNP calling


rna_star_variants <- function(
    fastq1,fastq2,


    ){

    source("~/iacsSVN/RScripts/pipelines/star/star_functions.R")
    tmp <- star2PassProc(fastq1=fastq1, fastq2=fastq2, outdir=outdir, STAR=star_exe, reffa=reffa,
                               genomeDir=reference_star,threads=cpu_star)
    cmds_star <- tmp$cmds
    rgid=sprintf("%s_%s", runid, lane, index) ## stays same even if demultplexing is wrong
    runPicard.fixRGTags(bam=bam, outBam=sorted_bam, javaMem=java_mem_rg, rgid=rgid,
                        , rgsm=sample)

}

if(FALSE){

    star_exe = "/scratch/rists/hpcapps/x86_64/STAR/2.3.0e/STAR"
    reference_star = "/scratch/rists/hpcapps/reference/human/broad_hg19/indexes/STAR/2.3.0e"
    cpu_star = 16
    reffa = "/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta"
    fastq1 = "fastq1"
    fastq2 = "fastq2"
    outdir = "."

}
