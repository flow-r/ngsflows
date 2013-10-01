#### PICARD, logFile sample place as out
## tool="EstimateLibraryComplexity.jar";params="R=/IACS1/NGS/hg19BWAIndex/Homo_sapiens_assembly19.fasta"
runPicard <- function(bamFile,tool,params="",outsuf,picardPath=getOption("ngs.picardPath"),outPath=getOption("ngs.sampleQAPath"),
                      java=getOption("ngs.java"),javaMem=getOption("ngs.javaMem"),insuf=".bam|.sam",force=FALSE,verbose=TRUE){
    outFile <- file.path(outPath,gsub(insuf,sprintf(".%s.out",outsuf),basename(bamFile)))
    logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    if(!force & file.exists(outFile))     return(outFile)
    cmd <- sprintf("time %s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/%s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT %s 2>>%s",
                   java,javaMem,picardPath,tool,bamFile,outFile,params,logFile)
    if(verbose) cat(cmd,"\n")
    system(cmd)
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}

runPicard.BamIndexStats <- function(bamFile,tool="BamIndexStats.jar",
                                    params="",outsuf,picardPath=getOption("ngs.picardPath"),outPath=getOption("ngs.sampleQAPath"),
                                    java=getOption("ngs.java"),javaMem=getOption("ngs.javaMem"),insuf=".bam|.sam",
                                    prefix=gsub(".bam|.sam","",basename(bamFile)),
                                    force=FALSE,
                                    verbose=TRUE){
  outFile <- file.path(outPath,gsub(insuf,sprintf(".%s.out",outsuf),basename(bamFile)))
  logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
  log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
  cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
  if(!force & file.exists(outFile))     return(outFile)
  cmd <- sprintf("time %s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/%s INPUT=%s VALIDATION_STRINGENCY=LENIENT %s 1>>%s 2>>%s",
                 java,javaMem,picardPath,tool,bamFile,params,outFile,logFile)
  if(verbose) cat(cmd,"\n")
  system(cmd)
  log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
  cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
  return(outFile)
}


### http://picard.sourceforge.net/command-line-overview.shtml picard is designed for 2gb ram
runPicard.fixRGTags <- function(bam,tool="AddOrReplaceReadGroups.jar",
                                picardPath=getOption("ngs.picardPath"),
                                outBam,rgid,rgsm,rglb,rgpu,rgpl="Illumina_HiSeq2000",rgcn='IACS',
                                javaMem=getOption("ngs.javaMem"),
                                javaTemp="/IACS1/tmp"){
    ## picard get the correct name
    cat("Working on ",bam,"\n")
    log=gsub(".bam",".log",outBam)
    cmd=sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/%s INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 1>> %s 2>> %s",
    getOption("ngs.java"),javaMem,javaTemp,picardPath,tool,bam,outBam,rgid,rglb,rgpl,rgpu,rgsm,rgcn,log,log)
    cat(cmd,"\n");system(cmd)
}


##debug(runSamFlagStat)
runSamFlagStat <- function(bamFile,samtools=file.path(getOption("ngs.samtoolPath"),"samtools"),outsuf="flagstat",
                           outPath=getOption("ngs.sampleQAPath"),insuf=".bam|.sam",verbose=TRUE,force=FALSE){
    dir.create(outPath,recursive=TRUE)
    outFile <- file.path(outPath,gsub(insuf,sprintf(".%s",outsuf),basename(bamFile)))
    logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN FLAGSTAT run\n")
    if(!force & file.exists(outFile))     return(outFile)
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    cmd <- paste("time ",samtools," flagstat ",bamFile," > ",outFile,sep="")
    if(verbose) cat(cmd,"\n")
    try(system(cmd))                  #spits out too much stuff...
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END FLAGSTAT run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}

## tool="DepthOfCoverage";params="-L MT";outsuf="dopMT";outPath="/IACS1/home/sseth/projects/moon/dopMT"
## tools="UnifiedGenotyper"; params="-L MT -glm BOTH -dfrac 1 --sample_ploidy 2 --dbsnp /IACS1/NGS/bundle1.2_b37/dbsnp_132.b37.vcf"
## type="gatk.SomaticIndel";tool <- "SomaticIndelDetector";logPath="/IACS1/GCC/log/indel";opts <- ""
runGatk <- function(bamFile,tool,params="",outsuf,gatkPath=getOption("ngs.gatkPath"),outPath=getOption("ngs.sampleQAPath"),
                    logPath=getOption("ngs.sampleQAPath"),
                    java=getOption("ngs.java"),javaMem=getOption("ngs.javaMem"),javaTemp=getOption("ngs.javaTemp"),
                    insuf=".bam|.sam",verbose=TRUE,force=FALSE){
    outFile <- file.path(outPath,gsub(insuf,sprintf(".%s",outsuf),basename(bamFile)))
    logFile <- file.path(logPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    if(!force & file.exists(outFile)){
        if(verbose) cat(sprintf("The outfile already exists: %s Use some force!",outFile));
        return(outFile)}
    cmd <- sprintf("time %s %s -Djava.io.tmpdir=%s -jar %s/GenomeAnalysisTK.jar -T %s -I %s -o %s %s 1>> %s 2>>%s",
                   java,javaMem,javaTemp,gatkPath,tool,bamFile,outFile,params,logFile,logFile)
    if(verbose) cat(cmd,"\n")
    try(system(cmd))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}

fastmean <- function(dat) {
    with(dat, sum(freq*value)/sum(freq) )
}
fastRMSE <- function(dat) {
    mu <- fastmean(dat)
    with(dat, sqrt(sum(freq*(value-mu)^2)/(sum(freq)-1) ) )
}
