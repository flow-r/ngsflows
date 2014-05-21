## includes annotations
## pairedmatfile: csv with samp, ref, sampbam, refbam
pipe_pindel <- function(pairedsamplefile,
                        flowname = "pindel",
                        samtools <- "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools",
                        pindel_exe = "/scratch/rists/hpcapps/x86_64/pindel/024t/pindel",
                        pindelpath = "/scratch/iacs/gcc/leveliii/pindel/tmp",
                        samtools_exe = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/samtools",
                        rscript_exe = "/scratch/rists/hpcapps/x86_64/R/3.1.0/bin/Rscript",
                        scripts = "/scratch/iacs/bin/pipelines/pipe.v2/scripts",
                        reffa ="/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta",
                        out_base_path="/scratch/iacs/gcc",
                        u = 0.1, execute=FALSE){
    mat <- read.csv(pairedsamplefile, as.is=TRUE)
    q_obj <- queue(type="torque",queue="iacs")
    for(i in 1:nrow(mat)){
        sample <- mat[i,"sample"]
        samplebam <- mat[i,"samplebam"]
        refbam <- mat[i,"refbam"]
        outprefix <- mat[i,"outprefix"]
        runid <- mat[i,"runid"]
        if(missing(pindelpath)) pindelpath = file.path(out_base_path,"leveliii",runid, "pindel") ## folder for samplemat
        configfile <- file.path("tmp", paste(outprefix, ".config", sep = ""))
        ## read_len <- system(sprintf("%s view %s | head -n 100 | awk 'BEGIN{}{len=len+length($10)}END{print len/100}'",
        ##                       samtools_exe, mat[i,"samplebam"]),intern=TRUE)
        read_len <- IACSUtil::getReadLength(samplebam, samtoolPath=dirname(samtools_exe))
        isize <- IACSUtil::getInsertSize(mat[i,"samplebam"], samtoolPath=dirname(samtools_exe))
        write.table(config, file = configfile, sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
        config <- rbind(c(samplebam, isize, gsub("_.*", "", basename(samplebam))),
                        c(refbam, isize, gsub("_.*", "", basename(refbam))))
        chrs <- system(sprintf("%s idxstats %s | cut -f1 | grep -v '*'",samtools_exe, samplebam),intern=TRUE)
        out_files <- sprintf("%s/%s_%s","tmp", outprefix, chrs)
        cmd_pindel <- sprintf("%s -f %s -i %s -u %s -T %s -o %s -c %s",
                       pindel_exe, reffa, configfile, u, cores, out_files, chrs)
        cmd_annotate_pindel <- sprintf("%s %s/runMyFunc.R annotatePindel tumorBAM_C=%s normalBAM_C=%s runID_C=NA pindel_C=%s_SI odir_C=%s oprefix_C=%s",rscript_exe, scripts,samplebam, refbam, out_files, "tmp", out_files)
        ## ----------------- Define the jobs
        j_obj_pindel <- job(q_obj=q_obj,cmds=cmd_aln,cpu=cpu_aln, submission_type="scatter", name="pindel",
                            next_job="pindel_annotate") ## no dependency
        j_obj_annotate <- job(q_obj=q_obj,cmds=cmd_annotate_pindel,cpu=cpu_aln, submission_type="scatter",
                              name="pindel_annotate", previous_job="pindel")


    }
}

if(FALSE){
    pairedsamplefile <- "/scratch/iacs/iacs_dep/sseth/tmp/pairedsamplefile.csv"
}


runPindel <- function(sampleBam, refBam, odir, oprefix,
    pindelPath = "/scratch/rists/hpcapps/x86_64/pindel/024t", gccPath = "/scratch/iacs/gcc",
    samtoolPath = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19", u = 0.1,
    commandOnly = FALSE, faName =
    "/scratch/rists/hpcapps/reference/human/hg19BWA0.7.5aIndex/Homo_sapiens_assembly19.fasta",
    cores = 4,force=FALSE){
    if(missing(odir)){
        odir <- file.path(gccPath, "leveliii", "pindel",
                             basename(dirname(sampleBam)))
    }
    if(missing(oprefix)){
    	bams <- ifelse(missing(sampleBam), NA, sampleBam)
    	bams <- c(bams, ifelse(missing(refBam), NA, refBam))
    	bams <- bams[!is.na(bams)]
       	oprefix <- paste(gsub("[\\._]+rg.*", "", basename(bams)), sep = "",
       	collapse = "__")
    }
    dir.create(odir, showWarning = FALSE, recursive = TRUE)
    configName <- file.path(odir, paste(oprefix, ".config", sep = ""))
    readLength <- getReadLength(sampleBam, samtoolPath = samtoolPath)
    sInsert <- ceiling(getInsertSize(sampleBam, samtoolPath = samtoolPath)) +
        2 * readLength
    outFile <- file.path(odir, oprefix)
    if(!missing(refBam)){
        readLength <- getReadLength(refBam, samtoolPath = samtoolPath)
        nInsert <- ceiling(getInsertSize(refBam, samtoolPath = samtoolPath)) +
            2 * readLength
        config <- rbind(c(sampleBam, sInsert, gsub("_.*", "", basename(sampleBam))),
                        c(refBam, nInsert, gsub("_.*", "", basename(refBam))))
        cmd <- paste(file.path(pindelPath, "pindel"), "-f ", faName, "-i",
                     configName,  "-u", u, "-T", cores, "-c ALL -o", outFile)
    }else{
        config <- t(c(sampleBam, sInsert, gsub("_.*", "",
                                               basename(sampleBam))))
        cmd <- paste(file.path(pindelPath, "pindel"), "-f ", faName, "-i", configName,
            "-T",  cores, "-u", u, "-c ALL -o", outFile)
    }
    write.table(config, file = configName, sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
    if(commandOnly){
    	return(cmd)
    }
    if(!file.exists(sprintf("%s_SI",outFile)) | force){ #file does not exists or force
        system(cmd)

    }else{
    	stop("File exits. Set force to TRUE if you want to overwrite the files")
    }
    return(sprintf("%s_SI",outFile))
}
