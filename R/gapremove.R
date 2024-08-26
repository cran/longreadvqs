#' Removing gap-rich positions and/or reads/sequences
#'
#' @description
#' Removes nucleotide positions (vertical) and/or reads/sequences (horizontal) that contain gaps more than the specified cut-off percentage from the alignment.
#'
#' @param fasta Input as a read or multiple sequence alignment in FASTA format
#' @param vgappct The percent cut-off of vertical gap (-), i.e., if a position in the alignment has %gap >= vgappct, that position will be removed.
#' @param hgappct The percent cut-off of horizontal gap (-), i.e., if a sequence or read in the alignment has %gap >= hgappct, that sequence or read will be removed.
#' @param fastaname Output file name in FASTA format
#'
#' @return FASTA read or multiple sequence alignment written out to the input directory
#' @export
#'
#' @examples
#' ## Locate input FASTA file-------------------------------------------------------------------------
#' fastafilepath <- system.file("extdata", "gaprichfast.fasta", package = "longreadvqs")
#'
#' ## Indicate output directory and file name---------------------------------------------------------
#' outfast <- tempfile()
#'
#' ## Remove positions with gap >= 60% and reads/sequences with gap >= 10%----------------------------
#' gapremove(fastafilepath, vgappct = 60, hgappct = 10, fastaname = outfast)
#'
#' @name gapremove

gapremove <- function(fasta, vgappct = 70, hgappct = 70, fastaname = "filteredfast.fasta"){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  seq <- readDNAStringSet(fasta)
  seq2 <- read.alignment(file = fasta, format = "fasta")

  if(missing(vgappct) & missing(hgappct)) {
    stop("Please specify vgappct or hgappct or both")
  }else if(missing(vgappct)){
    vgappct <- 100
  }else if(missing(hgappct)){
    hgappct <- 100
  }else{
    vgappct <- vgappct
    hgappct <- hgappct
  }

  #Vertical gap (position)
  frq <- seqinr::consensus(seq2, method = "profile")
  frq <- as.data.frame(frq)
  frqpc <- 100 * sweep(frq, 2, colSums(frq), `/`)
  frqpc <- as.data.frame(t(frqpc))
  hivgap <- apply(frqpc, 1, function(x) paste(names(which(x[1] >= vgappct)), collapse = ","))
  hivgap <- as.data.frame(hivgap)
  hivgap <- tibble::rownames_to_column(hivgap, "position")
  hivgap$position <- as.integer(hivgap$position)
  hivgap <- hivgap[!(is.na(hivgap$hivgap) | hivgap$hivgap==""), ]
  hivgap <- c(hivgap$position)
  orifast <- dss2df(seq)
  fastadf <- data.frame(str_split_fixed(orifast$seq, "", max(nchar(orifast$seq))))
  colnames(fastadf) <- c(1:ncol(fastadf))
  if(length(hivgap) == 0){
    newfast <- fastadf
  }else{
    newfast <- fastadf[,-which(names(fastadf) %in% hivgap)]
  }

  #Horizontal gap (sequence)
  hgap <- apply(X = newfast,MARGIN = 1,function(t){sum(grepl(pattern = "-",x = t,fixed = TRUE))})
  hgappc <- hgap*100/ncol(newfast)
  newfast2 <- newfast
  newfast2$hgappc <- hgappc
  rownames(newfast2) <- orifast$names
  if(length(which(newfast2$hgappc >= hgappct)) == 0){
    finalfastdf <- newfast2
  }else{
    finalfastdf <- newfast2[-(which(newfast2$hgappc >= hgappct)),]
  }
  finalfastdf <- finalfastdf[,-which(names(finalfastdf) %in% c("hgappc"))]
  finalfastdf <- data.frame(seq = do.call(paste, c(finalfastdf, sep="")), row.names = rownames(finalfastdf))
  finalfastdf <- tibble::rownames_to_column(finalfastdf, "seqname")

  writeFasta<-function(data, filename){
    fastaLines = c()
    for (rowNum in 1:nrow(data)){
      fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"seqname"], sep = "")))
      fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
    }
    fileConn<-file(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn)
  }
  writeFasta(finalfastdf, fastaname)
  message("#Total positions: ",nrow(frqpc),"\n","#Remaining positions: ",ncol(newfast),"\n","Removed position(s): ",paste(hivgap, collapse=", "),"\n",
      "#Total sequences: ",nrow(newfast),"\n","#Remaining sequences: ",nrow(finalfastdf),"\n","Removed sequence(s): ",paste(shQuote(row.names(newfast2)[which(newfast2$hgappc >= hgappct)], type="cmd"), collapse=", "))
}
