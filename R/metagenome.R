## metagenomics_flow

if(FALSE){
  fqs1 = ""
  flow_def = ""
}

meta_flow <- function(fqs1, fqs2, 
                      flow_def,
                      ### paths to folders
                      bwa_exe = "",
                      bwa_aln_opts = "",
                      bwa_sampe_opts="",
                      reflib = ""){
  sai_files1=file.path(gsub(".fastq.gz",".sai",basename(fqs1)))
  sai_files2=file.path( gsub(".fastq.gz",".sai",basename(fqs2)))
  bam_files=file.path(gsub(".fastq.gz",".bam",basename(fqs1)))
  cmd_aln1 <- sprintf("%s aln %s %s %s > %s",bwa_exe, bwa_aln_opts,reflib,fqs1, sai_files1)
  cmd_aln2 <- sprintf("%s aln %s %s %s > %s",bwa_exe, bwa_aln_opts,reflib,fqs2, sai_files2)
  cmd_sampe <- sprintf("%s sampe %s %s %s %s %s %s | %s view -Shu - > %s",
                       bwa_exe,bwa_sampe_opts,reflib,sai_files1,sai_files2,fqs1,fqs2,samtools_exe, bam_files)
  ## ----------------- START read group
  bamrg_files=sprintf("%s_rg.bam",bam_files)
  cmd_fixrg <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/AddOrReplaceReadGroups.jar INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT",
                       java_exe, java_mem, java_tmp, picard_dir, bam_files, bamrg_files, rgid, rglb, platform, rgpu, rgsm,
                       center)
  ## ----------------- START merging
  mergebam <- file.path(out_bampath,out_bam)
  bam_list <- paste("INPUT=", bamrg_files, sep = "", collapse = " ")
  java_mem <- "-Xmx8g";
  cmd_merge <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/MergeSamFiles.jar %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
                       java_exe, java_mem, java_tmp, picard_dir, bam_list, mergebam)
  return() 
}

