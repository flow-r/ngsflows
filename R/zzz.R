.onAttach <- function(lib, pkg){
  packageStartupMessage("ngsflows: genomic flows made faster")

  fls = unique(unlist(sapply(c("ngsflows"), flowr:::fetch_conf)))
  suppressMessages(flowr:::opts_flow$load(fls, check = FALSE))

  if(flowr:::opts_flow$get('verbose') > 1)
    packageStartupMessage("\nverbose level: 2 (debug mode)\nfollowing files are being loaded:\n", paste(fls, "\n"))

}
