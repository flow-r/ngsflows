if(FALSE){

  library(ngsflows)
  bam = "/rsrch2/iacs/ngs_runs/150107_SN1120_0324_AC6655ACXX/bams/DGomez-MOON2-Lung-127390-T_C6655ACXX-2-TGACCA.bwa_recalibed.bam"
  load_opts(fetch_conf("ngsflows.conf"), check = FALSE)
  out = bam_preprocess(bam, samplename = "127390-T", split_by_chr = TRUE)
  pip = fetch_pipes("bam_preprocess")
  
  set_opts(verbose = FALSE)
  def = as.flowdef(pip$def)
  
  fobj = to_flow(out$flowmat, pip$def)
  det = to_flowdet(fobj)
  
  fobj2 = submit_flow(fobj, execute = TRUE)
  
  
  plot_flow(def, pdffile = "inst/pipelines/bam_preprocess.pdf")
  ##testing on
  ## /scratch/iacs/iacs_dep/sseth/flows/AdenoidCysticSanger/PD3176a/merge-preproc-0dd3b13a-aef5-4cdc-bee6-15db06f2dc71/tmp
}
