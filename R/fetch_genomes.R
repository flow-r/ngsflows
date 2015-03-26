
#' @title genomes_fetch
#' @description genomes_fetch
#' @param base_url
#' @param species
#' @param src
#' @param build
#' @param ... passed onto \link{genomes_fetch}
#' @export
#' @importFrom RCurl url.exists
#' @examples \donotrun{
#' genomes_fetch(species = Homo_sapiens, src = 'NCBI', build = 'build37.2')
#' }
genomes_fetch <- function(genome_path = "~/flowr/genomes", ...){
  if(!file.exists(genome_path)) dir.create(genome_path, recursive = TRUE)
  setwd(genome_path)
  ## ---------- if any of these are missing call fetch to check
  tmp = genomes_avail(...)
  if(tmp$type == "tar"){
    url = tmp$url
    tarfl = tmp$lst
    message("Downloading tar...", tarfl)
    tmp <- getURL(url)
    message("Extracting tar...", tarfl)
    message("All one with", tarfl)
    invisible()
  }
}


#' @title genomes_avail
#' @description genomes_avail
#' @param species species
#' @param src src
#' @param build build
#' @param from from
#' @param base_url
#' @export
#' @importFrom RCurl getURL
#' @examples
#' gen = genomes_avail(species = 'Homo_Sapiens', from = 'igenomes')
genomes_avail <- function(species, src, build,
                          from = "igenomes",
                          verbose = FALSE,
                          base_url = "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com", ...){
  url = paste(base_url)
  head = "\n################################################\n"
  msg = c(head, 'Available Species:',head)
  type = "list"
  if(!missing(species)){
    msg = c(head, 'Available Sources:',head)
    url = paste(base_url, species, sep = "/")
    if(!missing(src)){
      msg = c(head, 'Available builds:',head)
      url = paste(base_url, species, src, sep = "/")
      if(!missing(build)){
        msg = c(head, 'Available files:',head)
        url = paste(base_url, species, src, build, sep = "/")
        type = "tar"
      }
    }
  }
  message(msg)
  if(verbose) message(url)
  lst = getURL(sprintf("%s/", url), dirlistonly = TRUE)
  message(lst)
  invisible(list(lst = lst, url = url, type = type))
}

