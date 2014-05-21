
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

    basecalls=sprintf("%s/intensity/%s/Data/Intensities/BaseCalls",Sys.getenv("GCC_PATH"),runid)
    basepath=sprintf("%s/intensity/%s",Sys.getenv("GCC_PATH"),runid)
    mask <- get_casava_mask(path=basepath)
    outpath=sprintf("%s/leveli/%s/plain",Sys.getenv("GCC_PATH"),runid)
    create_casava_plain_sheet(outPath=outpath, lanes=1:8, flowcellid=fcid)
    samplesheet=sprintf("%s/SampleSheetPlain.csv",outpath)

    casava="/scratch/rists/hpcapps/x86_64/CASAVA/1.8.2/bin/configureBclToFastq.pl"
    cmd=sprintf("%s --input-dir %s --output-dir %s --use-bases-mask %s --force --mismatches 1 --ignore-missing-stats --ignore-missing-bcl --ignore-missing-control --with-failed-reads --sample-sheet %s",casava,basecalls,output,mask,samplesheet)
    cat(cmd)
"cd $output;make -j 16 1>> casava.out 2>> casava.out"




}
