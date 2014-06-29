split.names.fastq <- function(files,format="$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"){
    ## process format:
    f <- gsub("\\$samplename\\$","(.*)",format)
    f <- gsub("\\$index\\$","([ATGC]*|NoIndex)",f)
    f <- gsub("\\$lane\\$","([0-9]*)",f)
    f <- gsub("\\$read\\$","([0-9]*)",f)
    f <- gsub("\\$num\\$","([0-9]*)",f)
    ## column names
    out=strsplit(format,"\\$")[[1]]
    cols=out[seq(2,length(out),by=2)]
    repl <- paste("\\",1:length(cols),sep="",collapse=",")
    mat <- gsub(f,repl,basename(files))
    mat <- data.frame(cbind(do.call(rbind,strsplit(mat,",")),files))
    colnames(mat) <- c(cols,"files")
    return(mat)
}
split_names_fastq=split.names.fastq

create_sample_mat <- function(path, project, subproject, runid, outpath, format,
                              pattern = "fastq.gz|fq.gz|fastq|fq"){
    fqs <- list.files(path,pattern = pattern,full.names=TRUE,recursive=TRUE)
    if(missing(project)) project = basename(path)
    if(missing(subproject)) subproject = substr(project, 1, 2)
    if(missing(outpath)) outpath = "." ## folder for samplemat
    if(missing(runid)) runid = basename(dirname(dirname(dirname(fqs[1])))) ## runid
    if(missing(format)){
        if(grepl("_S1.*fastq.gz",fqs[1])){ ## miseq output
            ## ------ casava output
            format <- "$samplename$_S[0-9]*_L00$lane$_R$read$_$num$.fastq.gz"
        }else if(grepl(".*_([ATGC]*|NoIndex).*L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz",basename(fqs[1]))){
            ## ------ casava output
            format <- "$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"
        }
    }
    fq_mat <- split_names_fastq(files = fqs, format = format)
    ## ------- cleanup things
    fq_mat$samplename <- gsub("_","",fq_mat$samplename)
    out_basename <- sprintf("%s-%s-%s_%s", project, subproject, fq_mat$samplename, runid)
    sorted_bam <- sprintf("%s_rg.sorted.bam",out_basename)
    recal_bam <- sprintf("%s_rg.sorted.recalibed.bam", out_basename)
    fq_mat <- cbind(fq_mat, out_basename, sorted_bam, recal_bam, runid, project, subproject)
    if(!file.exists(outpath)) dir.create(outpath)
    write.csv(fq_mat, sprintf("%s/%s_%s_%s_sample_mat.csv", outpath, project, subproject,
                              runid), row.names=FALSE)
    return(fq_mat)
}

check_fastq_sheet <- function(mat){
    cat("There are", length(unique(mat$samplename)), "samples in these dataset\n")
    dat_list <- split.data.frame(mat, mat$samplename)
    tmp <- sapply(1:length(dat_list), function(i){
        tmp <- tapply(dat_list[[i]]$files, dat_list[[i]]$read, length)
        if(diff(tmp) != 0) stop("Number of fastq files are not the same in",
                   dat_list[[i]]$samplename[1], "sample")
        return(0)
    })
}

read_sample_sheet <- function(x){
    ext <- file_ext(x)
    if(ext=="tsv"){
        mat <- read.table(x, as.is=TRUE, sep="\t", header=TRUE)
    }else if(ext=="csv"){
        mat <- read.csv2(x, as.is=TRUE, comment.char = '#', strip.white=TRUE,
                         blank.lines.skip=TRUE, sep=",", header=TRUE)
    }
    check_fastq_sheet(mat)
    return(mat)
}

if(FALSE){
    path = "/scratch/iacs/iacs_dep/sseth/data/mt_scc/Project_VC_SCC-NI56-1"
    runid= "140626_M01692_0046_000000000-A99M8"
    runid="140627_M01692_0047_000000000-A9J0V";
    #debug(create_sample_mat)

    path=sprintf("/scratch/iacs/gcc/leveli/%s",runid)
    fq_mat <- create_sample_mat(path,project="ANDY",subproject="Futreal-AS", runid=runid,
                                outpath="~/flows/ANDY-Futreal-AS")


}
