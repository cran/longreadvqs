#' Filtering highly dissimilar reads/sequences out of the alignment
#'
#' @description
#' Removes reads/sequences of which Hamming similarity to the consensus of all reads/sequences in the alignment is less than the specified quantile (qt) of the similarity distribution.
#'
#' @param fasta Input as a read or multiple sequence alignment in FASTA format
#' @param qt If Hamming similarity score of a read/sequence to the consensus of all reads/sequences is less than the specified quantile (qt) of the similarity distribution, that read/sequence will be removed.
#' @param fastaname Output file name in FASTA format
#'
#' @return FASTA read or multiple sequence alignment written out to the input directory
#' @export
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings width
#' @importFrom seqinr read.alignment
#' @importFrom seqinr consensus
#' @importFrom stringdist stringsim
#' @importFrom stats quantile
#'
#' @examples
#' ## Locate input FASTA file-------------------------------------------------------------------------
#' fastafilepath <- system.file("extdata", "dissimfast.fasta", package = "longreadvqs")
#'
#' ## Indicate output directory and file name---------------------------------------------------------
#' outfast <- tempfile()
#'
#' ## Remove reads/sequences that the similarity < 1st quartile (0.25 quantile)-----------------------
#' filtfast(fastafilepath, qt = 0.25, fastaname = outfast)
#'

filtfast <- function(fasta, qt = 0.25, fastaname = "filteredfast.fasta"){
  dss2df <- function(dss) data.frame(width= Biostrings::width(dss), seq=as.character(dss), names=names(dss))
  writeFasta <- function(data, filename){
    fastaLines = c()
    for (rowNum in 1:nrow(data)){
      fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"names"], sep = "")))
      fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
    }
    fileConn<-file(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn)
  }
  seq <- readDNAStringSet(fasta)
  seq2 <- read.alignment(file = fasta, format = "fasta")
  fastadf <- dss2df(seq)
  consensus <- seqinr::consensus(seq2, method = "majority")
  consensus <- paste(consensus, collapse = '')
  fastadf$sim <- sapply(fastadf$seq, function(x){stringdist::stringsim(x, consensus, method = "hamming")})
  fastafilt <- fastadf[fastadf$sim >= quantile(fastadf$sim, qt), ]
  writeFasta(fastafilt, fastaname)
  message("#Total sequences: ",nrow(fastadf),"\n","#Remaining sequences: ",nrow(fastafilt),"\n", "Removed sequence(s): ",paste(shQuote(row.names(fastadf)[which(fastadf$sim < quantile(fastadf$sim, qt))], type="cmd"), collapse=", "))
}
