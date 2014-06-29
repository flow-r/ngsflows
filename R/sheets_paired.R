
require(tools)
read_paired_samplesheet <- function(x, ...){
    ext <- file_ext(x)
    if(ext=="tsv"){
        mat <- read.table(x, as.is=TRUE, sep="\t", header=TRUE)
    }else if(ext=="csv"){
        mat <- read.csv2(x, as.is=TRUE, comment.char = '#', strip.white=TRUE,
                         blank.lines.skip=TRUE, sep=",", header=TRUE)
    }
    mat <- mat[! (is.na(mat$samplename) | mat$samplename %in% c("NA", "NULL", "")), ]
    ## convert paired_mat for xenome, calling snps etc
    single_mat <- rbind(cbind(mat$project, mat$samplename, mat$sampbam, mat$db_sampleid),
                        cbind(mat$project, mat$refname, mat$refbam, mat$db_refid))
    colnames(single_mat) <- c("project", "samplename", "bam", "db_id")
    single_mat <- data.frame(unique(single_mat),stringsAsFactors=FALSE)
    return(list(single_mat=single_mat, paired_mat=mat))
}

tooling_paired_samplesheet <- function(x, outfile, tumor.only=FALSE, normal.only=FALSE){
    colnames(x)  <- tolower(colnames(x))
    project=apply(x[,c("project","subproject")], 1, paste, collapse="-")
    if(tumor.only){
        out_mat <- data.frame(project, samplename=x$tumorsampleid, refname=NA,
                              sampbam=x$tumorsamplepreprocbam, refbam=NA,
                              db_sampleid=x$tumorsamplerunitemid, db_refid=0)
    }
    if(!missing(outfile)){
        dir.create(dirname(outfile),recursive=TRUE)
        write.table(out_mat, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    }
    return(out_mat)
}
