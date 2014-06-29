
if(FALSE){

    require(ibm.ppm.tools)
    x <- "/scratch/iacs/iacs_dep/sseth/flows/GIULIO-GSC-JG/paired_samplesheet.csv"
    ##apply(mat, 1, function(z) getPairOutPrefixBam(tbam=z["sampbam"], nbam=z["refbam"]))
    .submit_fastq_preprocess(x = x, species = "human",
                    outpath = "/scratch/iacs/gcc/leveliii/mutect/130906_SN1120_0275_AC2A7CACXX")

}

.submit_fastq_preprocess <- function(x, species, outpath){
    ## source("~/Dropbox/public/github_ngsutils/R/sample_sheets.R")
    ## debug(read_paired_samplesheet)
    mat <- read_sample_sheet(x)

    for(i in 1:nrow(mat)){
        out_basename <- sprintf("%s-%s-%s_%s",s.project,s.subproject,sample,s.runid)
        sorted_bam <- sprintf("%s/%s_rg.sorted.bam", bampath, out_basename)
        recal_bam <- sprintf("%s/%s_rg.sorted.recalibed.bam", bampath, out_basename)
        fqs1=as.c(subset(sample.mat,read==1)$files)
        fqs2=as.c(subset(sample.mat,read==2)$files)
        o.sorted <- submit.fastq.sorted.bam(fqs1=fqs1,fqs2=fqs2,outBam=out_basename,
                                            out.bampath=bampath,sample=sample,se.or.pe=se.or.pe,
                                            bwa.aln.opts="-l 50 -n 6 -k 4",runid=runid,
                                            lane=s.lane,platform='illumina',center='IACS',
                                            species=species, addl.param=addl.param)
        ##o.sorted$flowid=0
        o.recal <- submit.sorted.recal.bam(inbam= sorted.bam,out.bampath=bamPath,species=species,
                                           trigger.flow.ids=o.sorted$flowid,
                                           cleanup=0)
    }
}


