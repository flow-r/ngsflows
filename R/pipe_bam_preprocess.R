bam_preprocess <- function(inbam, outbam,
                               java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                               java_mem = "-Xmx8g",
                               java_tmp = "tmp",
                               gatk_jar = "/scratch/rists/hpcapps/x86_64/gatk/3.1-1/GenomeAnalysisTK.jar",
                               picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112",
                               reffa = "/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta",
                               gatk_realigner_opt = "--known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/mills_and_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf --known /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores",
                               gatk_baserecalibrator_opt  = "-knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/dbsnp/dbsnp_138.b37.excluding_sites_after_129.vcf -knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/mills_and_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites /scratch/iacs/iacs_dep/sseth/reference/human/b37/annotations/1000g/1000G_phase1.indels.b37.vcf --allow_potentially_misencoded_quality_scores -nct 16"){
    dedupbam <- gsub(".bam",".marked.bam",basename(inbam))
    metricsfile <- gsub(".bam",".marked.metrics",basename(inbam))
    cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/MarkDuplicates.jar INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORTED=true   VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                           java_exe,java_mem,java_tmp,picard_dir,inbam,dedupbam,metricsfile)
    ## ------------ realign
    intervalsfile <- gsub(".bam",".realign.tab",basename(inbam))
    realignedbam <- gsub(".bam",".realigned.bam",basename(inbam))
    cmd_targetcreater <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator -R %s -I %s -o %s%s",
                                 java_exe,java_mem,java_tmp, gatk_jar,reffa,dedupbam,intervalsfile, gatk_realigner_opt)
    cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s",
                           java_exe,java_mem, gatk_jar,reffa,inbam,intervalsfile,realignedbam, gatk_realigner_opt)
    ## ------------ base recalibration
    recalibbam <- gsub(".bam",".recalibed.bam",basename(inbam))
    recalibtabfile <- gsub(".bam",".recalib.tab",basename(inbam))
    cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -R %s -I %s -o %s %s",
                               java_exe,java_mem,java_tmp, gatk_jar,reffa,realignedbam,recalibtabfile, gatk_baserecalibrator_opt)
    cmd_printreads <- sprintf("%s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/GenomeAnalysisTK.jar -T PrintReads -R %s -I %s -BQSR %s -o %s -allowPotentiallyMisencodedQuals",
                              java_exe,java_mem, gatk_jar,reffa,realignedbam,recalibtabfile,recalibbam)
    cmds <- c(cmd_markdup=cmd_markdup, cmd_targetcreater=cmd_targetcreater, cmd_realign=cmd_realign,
              cmd_baserecalib=cmd_baserecalib, cmd_printreads=cmd_printreads)
    return(cmds)
}
