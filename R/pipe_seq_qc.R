## A class that contains bam file information
## Copyright 2014, Sahil Seth, all rights reserved
## sahil.seth@me.com
## A few functions to supplement those already in this package.
#### -----------------------

## ---------- run this to set all params:
if(FALSE){
    runid="140506_SN1120_0305_AC3ERWACXX"
    fcid="AC3ERWACXX"
}

pipe_qc <- function(runid, check_barcode=TRUE, aln=TRUE, rnaseq_qc=FALSE, capture_qc=FALSE, capture_qc_report=FALSE){
    ## takes different routes for miseq and hiseq
    ## this would do fastqc
    ## if check_barcode is true make a plain run and get list of IDs
    ## if aln is true: run alignment using bwa
    ## if rnaseq_qc is true; run that
    #system("source /scratch/iacs/bin/pipe.v2/scripts/set.environ")
    ## ------------ CASAVA
    basecalls=sprintf("%s/intensity/%s/Data/Intensities/BaseCalls",Sys.getenv("GCC_PATH"),runid)
    basepath=sprintf("%s/intensity/%s",Sys.getenv("GCC_PATH"),runid)
    mask <- get_casava_mask(path=basepath)
    outpath=sprintf("%s/leveli/%s/plain",Sys.getenv("GCC_PATH"),runid)
    create_casava_plain_sheet(outPath=outpath, lanes=1:8, flowcellid=fcid)
    samplesheet=sprintf("%s/SampleSheetPlain.csv",outpath)
    casava_opts <- "--mismatches 1 --ignore-missing-stats --ignore-missing-bcl --ignore-missing-control --with-failed-reads"
    casava="/scratch/rists/hpcapps/x86_64/CASAVA/1.8.2/bin/configureBclToFastq.pl"
    cmd_casava <- sprintf("%s --input-dir %s --output-dir %s --use-bases-mask %s --sample-sheet %s --force %s",
        casava,basecalls,output,mask,samplesheet, casava_opts)
    cat(cmd_casava)
    cmd_casava_make <- sprtintf("cd %s;make -j 16 1>> casava.out 2>> casava.out",outpath)
    cat(cmd_casava_make)
    count_indexes="/scratch/iacs/iacs_dep/sseth/Dropbox/public/github.r-ngs-utils/inst/files/count_fastq_indexes.sh"
    cmd_getindextable <- sprintf("cores=%s %s %s/Project_SampleProject/Sample_Sample%s", 16, count_indexes,fqpath,lanes)
    cmd_email_report <- sprintf("cd %s;head Project_SampleProject/Sample_Sample*/total_summary.txt | sendit \"%s\"",
                                fqpath,fqpath)
    system(cmd_email_report)

}
