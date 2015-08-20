mypack:::reload('flow')

pipe_bismark <- function(
  fq_path = "/scratch/iacs/iacs_dep/sseth/data/giulio_pushan/Project_ME_PD-1_PM38/Sample_6-27_GFP",
  out_path = "/rsrch1/iacs/iacs_dep/sseth/data"
  perl_exe = "/risapps/rhel6/perl/5.10.1/bin/perl",
  bismark_dir = "/scratch/iacs/iacs_dep/sseth/apps/bismark_v0.12.5",
  bismark_index_path = "/scratch/iacs/data/reference/human/broad_hg19/indexes/bismark",
  bowtie_path = "/risapps/rhel5/bowtie/1.0.1",
  samtools_exe = "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/bin/samtools",
  java_exe = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
  java_mem = "-Xmx4g",
  picard_path = "/scratch/rists/hpcapps/x86_64/picard/1.112",
  rscript_exe = "/scratch/rists/hpcapps/rhel6/R/3.1.0/bin/Rscript",
  trim_galore_exe = "/scratch/iacs/iacs_dep/sseth/apps/trim_galore/trim_galore",
  q_type = "lsf",
  q_queue = "normal",
  flowname = "bismark_pipe",
  flow_base_path = "/scratch/iacs/iacs_dep/sseth/flows",
  project = "giulio_pushan",
  execute = TRUE, submit = TRUE,
  cpu_trim_galore = 1, cpu_bismark = 4, cpu_merge_bam = 2, cpu_meth_extractor = 2, cpu_sort = 1, cpu_methylseq = 1){
  fq_files = list.files(fq_path, pattern = "fastq.gz$", full.names = TRUE)
  cmds.trim_galore = sprintf("%s --quality 20 --phred33 --fastqc --rrbs --output_dir ./ --fastqc_args '--outdir ./ --casava' %s",
                             trim_galore_exe, fq_files)
  qual_fq_files = gsub(".fastq.gz", "_trimmed.fq.gz", basename(fq_files)); #file.exists(qual_fq_files)
  cmds.bismark = sprintf("%s %s/bismark -n 1 -l 36 --path_to_bowtie %s --samtools_path %s -bam %s %s", 
                         perl_exe, bismark_dir, bowtie_path, samtools_exe, bismark_index_path, qual_fq_files)
  bam_list = gsub("fastq.gz$", "fastq.gz_bismark.bam", basename(fq_files))
  sample = gsub("(.*)(Sample_)(.*)/(.*)", "\\3", fq_files[1])
  merge_bam = basename(gsub("(.*)(Sample_)(.*)/(.*)", "\\1\\2\\3/\\3_bismark_merged.bam", fq_files[1]))
  assume_sorted = FALSE
  cmd.merge = sprintf("%s %s -jar %s/MergeSamFiles.jar %s OUTPUT=%s ASSUME_SORTED=%s VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false",
                      java_exe, java_mem, picard_path, paste("INPUT=", bam_list, sep = "", collapse = " "),
                      merge_bam, assume_sorted)
  sorted_prefix = gsub("bam$", "sorted", merge_bam)
  sorted_sam = gsub("bam$", "sorted.sam", merge_bam)
  cmd.sort = sprintf("%s sort %s %s; %s view %s.bam > %s", 
                     samtools_exe, merge_bam, sorted_prefix,
                     samtools_exe, sorted_prefix, sorted_sam)
  cmd.methylseq = sprintf("%s /scratch/iacs/iacs_dep/sseth/bin/runMyFunc.R read.bismark location_C=%s sample.id_C=%s assembly_C=hg19 save.folder_C=. save.context_C=CpG,CHG,CHH",
                          rscript_exe, sorted_sam, sample)
  cmd.meth_extractor = sprintf("%s/bismark_methylation_extractor -s --comprehensive %s",
                               bismark_dir, merge_bam)
  cmd.transfer = sprintf("transfer_files.sh -o %s %s *_CpG.txt",
                         out_path, merge_bam)
  require(flow)
  q_obj <- queue(type=q_type,queue=q_queue)
  j_obj_trim_galore <- job(q_obj = q_obj, name = "trim_galore", cmds = cmds.trim_galore, cpu = cpu_trim_galore, 
                           submission_type="scatter")
  j_obj_bismark <- job(q_obj = q_obj, name = "bismark", cmds = cmds.bismark, cpu = cpu_bismark, submission_type="scatter",
                       previous_job = "trim_galore", dependency_type = "serial")
  j_obj_merge_bam <- job(q_obj = q_obj, name = "merge_bam", cmds = cmd.merge, cpu = cpu_merge_bam,
                         submission_type = "serial", previous_job = "bismark", dependency_type = "gather")
  j_obj_meth_extractor <- job(q_obj = q_obj, name = "meth_extractor", cmds = cmd.meth_extractor, cpu = cpu_meth_extractor,
                              submission_type = "serial", previous_job = "merge_bam", dependency_type = "serial")
  j_obj_sort <- job(q_obj = q_obj, name = "sort", cmds = cmd.sort, cpu = cpu_sort,
                              submission_type = "serial", previous_job = "merge_bam", dependency_type = "serial")
  j_obj_methylseq <- job(q_obj = q_obj, name = "methylseq", cmds = cmd.methylseq, cpu = cpu_methylseq,
                    submission_type = "serial", previous_job = "sort", dependency_type = "serial")
  j_obj_transfer <- job(q_obj = queue(q_obj,queue = 'filetransfer'), name = "transfer", cmds = cmd.transfer, 
                        cpu = 1, submission_type = "serial", previous_job = "meth_extractor",
                        dependency_type = "serial")
  f_obj <- flow(jobs = list(j_obj_trim_galore, j_obj_bismark, j_obj_merge_bam, j_obj_meth_extractor, j_obj_sort, 
                            j_obj_methylseq, j_obj_transfer),
                name = flowname, mode = "scheduler", flow_base_path = flow_base_path, 
                desc = sprintf("%s/%s-%s", project, flowname, sample))
  #plot_flow(f_obj)
  if(submit) submit_flow(f_obj, execute = execute)
  return(f_obj)
}