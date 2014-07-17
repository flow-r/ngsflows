
processParams <- function(file,pkgname="Rngsutils",verbose=FALSE){
  params <- read.table(file,as.is=TRUE,sep="\t",strip.white=TRUE)
  for(i in 1:dim(params)[1]){
    l=list(params[i,2])
    names(l)=params[i,1]
    if(grepl("file.ngs",names(l))){
      l <- file.path(path.package(pkgname),"files",paste(params[i,2]))
      names(l) <- substr(paste(params[i,1]),6,nchar(params[i,1]))
      l <- as.list(l)
    }
    options(l)
  }
  if(verbose)
    print(params)
  return(params)
}

.onAttach <- function(libname, pkgname){
  params <- try(processParams(file.path(libname, pkgname, "files","default.params"),pkgname,verbose=FALSE))
}
