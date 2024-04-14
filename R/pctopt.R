#' Optimizing cut-off percentage for error minimization
#'
#' @description
#' Finds an optimal cut-off percentage for error minimization (in vqssub, vqsassess, and vqscustompct functions) that can decrease the number of singleton haplotypes to less than the desired percentage of the total reads.
#'
#' @param fasta Input as a read alignment in FASTA format
#' @param pctsing The desired percentage of singleton haplotypes relative to the total reads in the alignment.
#' @param method Sequencing error minimization methods that replace low frequency nucleotide base (less than the "pct" cut-off) with consensus base of that position ("conbase": default) or with base of the dominant haplotype ("domhapbase").
#' @param samplingfirst Downsampling before (TRUE) or after (FALSE: default) the error minimization.
#' @param gappct The percent cut-off particularly specified for gap (-). If it is not specified or less than "pct", "gappct" will be equal to "pct" (default).
#' @param ignoregappositions Replace all nucleotides in the positions in the alignment containing gap(s) with gap. This will make such positions no longer single nucleotide variant (SNV). The default is "FALSE".
#' @param samsize Sample size (number of reads) after down-sampling. If it is not specified or more than number of reads in the original alignment, down-sampling will not be performed (default).
#' @param label String within quotation marks indicating name of read alignment (optional).
#'
#' @return An optimal cut-off percentage for error minimization of an input sample and parameter settings. If label is specified, the output will be a data frame with percentage of singleton haplotypes at each cut-off percentage from zero to the optimal cut-off percentage.
#' @export
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings width
#' @importFrom Biostrings nmismatch
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom seqinr read.alignment
#' @importFrom seqinr as.alignment
#' @importFrom seqinr consensus
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import QSutils
#'
#' @examples
#' ## Locate input FASTA file-------------------------------------------------------------------------
#' fastafilepath <- system.file("extdata", "s1.fasta", package = "longreadvqs")
#'
#' ## Find an cut-off percentage that creates singleton haplotypes less than 50% of the alignment.----
#' pctopt(fastafilepath, pctsing = 50, label = "s1")
#'
#' @name pctopt

utils::globalVariables("newcol")

pctopt <- function(fasta, pctsing = 0, method= c("conbase", "domhapbase"), samplingfirst = FALSE, gappct = 50, ignoregappositions = FALSE, samsize = 100, label = "sample"){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  seq <- readDNAStringSet(fasta)
  seq2 <- read.alignment(file = fasta, format = "fasta")
  fulldepth <- length(seq)

  if(samsize >= length(seq)){
    samplingwhen <- "no_sampling"
  }else if(samplingfirst == TRUE){
    samplingwhen <- "before"
  }else if(samplingfirst == FALSE | missing(samplingfirst)) {
    samplingwhen <- "after"
  }else{
    samplingwhen <- "after"
  }

  method <- match.arg(method)

  if(method == "conbase"){
    ##sub-sampling before
    if(samplingfirst == FALSE | missing(samplingfirst) | missing(samsize)){
      seq <- seq
      seq2 <- seq2
    }else if(samplingfirst == TRUE){
      if(missing(samsize)){
        seq <- seq
        seq2 <- seq2
      }else if(samsize < length(seq)){
        samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
        seq <- seq[names(seq) %in% samreads]
        seq3 <- seqinr::as.alignment()
        seq3$nb <- samsize
        seq3$nam <- seq2$nam[seq2$nam %in% samreads]
        seq3$seq <- seq2$seq[seq2$nam %in% samreads]
        seq2 <- seq3
      }else{
        seq <- seq
        seq2 <- seq2
      }
    }else{
      seq <- seq
      seq2 <- seq2
    }

    frq <- seqinr::consensus(seq2, method = "profile")
    frq <- as.data.frame(frq)
    frqpc <- 100 * sweep(frq, 2, colSums(frq), `/`)
    frqpc <- as.data.frame(t(frqpc))

    store <- data.frame()

    for(pct in 0:100) {
      lowfrq <- apply(frqpc, 1, function(x) paste(names(which(x[1:5] < pct & x[1:5] > 0)), collapse = ","))
      lowfrq <- as.data.frame(lowfrq)

      if(gappct > pct){
        lowfrq2 <- apply(frqpc, 1, function(x) paste(names(which(x[1] < gappct & x[1] > pct)), collapse = ","))
        lowfrq2 <- as.data.frame(lowfrq2)
        lowfrq3 <- cbind(lowfrq, lowfrq2)
        lowfrq4 <- lowfrq3 %>% mutate(newcol = case_when(lowfrq2 == "" & lowfrq != "" ~ lowfrq, lowfrq2 != "" & lowfrq != "" ~ paste(lowfrq2, lowfrq, sep = ","), lowfrq2 != "" & lowfrq == "" ~ lowfrq2)) %>% select(newcol)
        colnames(lowfrq4) <- "lowfrq"
        lowfrq <- lowfrq4
      }else if(missing(gappct)) {
        lowfrq <- lowfrq
      }else{
        lowfrq <- lowfrq
      }

      lowfrq <- tibble::rownames_to_column(lowfrq, "position")
      lowfrq$position <- as.integer(lowfrq$position)
      lowfrq <- separate_rows(lowfrq,lowfrq,sep=",")
      lowfrq <- as.data.frame(lowfrq)
      lowfrq <- lowfrq[!(is.na(lowfrq$lowfrq) | lowfrq$lowfrq==""), ]
      lowfrq$lowfrq = toupper(lowfrq$lowfrq)
      maxfrq <- apply(frqpc, 1, function(x) paste(names(which(x==max(x)))))
      maxfrq <- lapply(maxfrq, function(x) replace(x, length(x) != 1, "-"))
      maxfrq <- sapply(maxfrq,"[[",1)
      maxfrq <- as.data.frame(maxfrq)
      maxfrq <- tibble::rownames_to_column(maxfrq, "position")
      maxfrq$position <- as.integer(maxfrq$position)
      maxfrq$maxfrq = toupper(maxfrq$maxfrq)
      if(dim(lowfrq)[1] == 0){
        message("No low frequency SNV")
      }else{
        lowfrq <- merge(x = lowfrq, y =  maxfrq, by = "position", all.x = TRUE)
      }

      dfseq <- dss2df(seq)
      dfseq <- data.frame(str_split_fixed(dfseq$seq, "", max(nchar(dfseq$seq))))

      for (i in 1:nrow(lowfrq)){
        if(dim(lowfrq)[1] == 0){
          message("No low frequency SNV")
        }else{
          dfseq[,lowfrq[i,1]][dfseq[,lowfrq[i,1]] == lowfrq[i,2]] <- lowfrq[i,3]
        }
      }

      if(ignoregappositions == TRUE){
        gapfound <- apply(dfseq, 1, function(x) paste(names(which(x=="-"))))
        gapposition <- unique(unlist(gapfound))
        for (i in gapposition){
          dfseq[[i]] <- "-"
        }
      }else if(ignoregappositions == FALSE | missing(ignoregappositions)){
        dfseq <- dfseq
      }else{
        dfseq <- dfseq
      }

      dfseq2 <- unite(dfseq, col='seq', c(names(dfseq[1:ncol(dfseq)])), sep='')
      readnames <- dss2df(seq)
      dfseq2$names <- readnames$names
      seq <- DNAStringSet(c(dfseq2$seq))
      seq@ranges@NAMES <- c(dfseq2$names)

      ##sub-sampling after
      if(samplingfirst == TRUE | missing(samsize)){
        seq <- seq
      }else if(samplingfirst == FALSE | missing(samplingfirst)){
        if(samsize < length(seq)){
          samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
          seq <- seq[names(seq) %in% samreads]
        }else{
          seq <- seq
        }
      }else{
        if(samsize < length(seq)){
          samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
          seq <- seq[names(seq) %in% samreads]
        }else{
          seq <- seq
        }
      }

      ##vqs analysis
      hap <- QSutils::Collapse(seq)
      if(length(hap$hseqs) <= 2){
        hapre <- hap
        names(hapre)[2] <- "seqs"
      }else{
        hapcor <- CorrectGapsAndNs(hap$hseqs[2:length(hap$hseqs)],hap$hseqs[[1]])
        hapcor <- c(hap$hseqs[1],hapcor)
        hapre <- QSutils::Recollapse(hapcor,hap$nr)
      }

      depth <- length(seq)
      nsingleton <- sum(hapre$nr == 1)
      pctsingleton <- nsingleton*100/depth

      storeout <- data.frame(label, pct, pctsingleton)
      store <- rbind(store, storeout)

      if(pctsingleton <= pctsing) {
        break
      }
    }
    if(missing(label)){
      return(pct)
    }else{
      return(store)
    }
  }

  if(method == "domhapbase"){
    ##sub-sampling before
    if(samplingfirst == FALSE | missing(samplingfirst) | missing(samsize)){
      seq <- seq
      seq2 <- seq2
    }else if(samplingfirst == TRUE){
      if(missing(samsize)){
        seq <- seq
        seq2 <- seq2
      }else if(samsize < length(seq)){
        samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
        seq <- seq[names(seq) %in% samreads]
        seq3 <- seqinr::as.alignment()
        seq3$nb <- samsize
        seq3$nam <- seq2$nam[seq2$nam %in% samreads]
        seq3$seq <- seq2$seq[seq2$nam %in% samreads]
        seq2 <- seq3
      }else{
        seq <- seq
        seq2 <- seq2
      }
    }else{
      seq <- seq
      seq2 <- seq2
    }

    frq <- seqinr::consensus(seq2, method = "profile")
    frq <- as.data.frame(frq)
    frqpc <- 100 * sweep(frq, 2, colSums(frq), `/`)
    frqpc <- as.data.frame(t(frqpc))

    store <- data.frame()

    for(pct in 0:100) {
      lowfrq <- apply(frqpc, 1, function(x) paste(names(which(x[1:5] < pct & x[1:5] > 0)), collapse = ","))
      lowfrq <- as.data.frame(lowfrq)

      if(gappct > pct){
        lowfrq2 <- apply(frqpc, 1, function(x) paste(names(which(x[1] < gappct & x[1] > pct)), collapse = ","))
        lowfrq2 <- as.data.frame(lowfrq2)
        lowfrq3 <- cbind(lowfrq, lowfrq2)
        lowfrq4 <- lowfrq3 %>% mutate(newcol = case_when(lowfrq2 == "" & lowfrq != "" ~ lowfrq, lowfrq2 != "" & lowfrq != "" ~ paste(lowfrq2, lowfrq, sep = ","), lowfrq2 != "" & lowfrq == "" ~ lowfrq2)) %>% select(newcol)
        colnames(lowfrq4) <- "lowfrq"
        lowfrq <- lowfrq4
      }else if(missing(gappct)) {
        lowfrq <- lowfrq
      }else{
        lowfrq <- lowfrq
      }

      lowfrq <- tibble::rownames_to_column(lowfrq, "position")
      lowfrq$position <- as.integer(lowfrq$position)
      lowfrq <- separate_rows(lowfrq,lowfrq,sep=",")
      lowfrq <- as.data.frame(lowfrq)
      lowfrq <- lowfrq[!(is.na(lowfrq$lowfrq) | lowfrq$lowfrq==""), ]
      lowfrq$lowfrq = toupper(lowfrq$lowfrq)
      domhap <- QSutils::Collapse(seq)
      if(domhap$nr[1] == 1){
        stop("No dominant haplotype, please use conbase method instead")
      }else{
        domhap <- dss2df(domhap$hseqs)
      }

      domhap <- domhap["seq"]
      domhap <- data.frame(str_split_fixed(domhap$seq, "", max(nchar(domhap$seq))))
      domhap <- domhap[1,]
      domhap <- as.data.frame(t(domhap))
      domhap$position <- 1:nrow(domhap)
      colnames(domhap) <- c("domhap", "position")

      if(dim(lowfrq)[1] == 0){
        message("No low frequency SNV")
      }else{
        lowfrq <- merge(x = lowfrq, y =  domhap, by = "position", all.x = TRUE)
      }

      dfseq <- dss2df(seq)
      dfseq <- data.frame(str_split_fixed(dfseq$seq, "", max(nchar(dfseq$seq))))

      for (i in 1:nrow(lowfrq)){
        if(dim(lowfrq)[1] == 0){
          message("No low frequency SNV")
        }else{
          dfseq[,lowfrq[i,1]][dfseq[,lowfrq[i,1]] == lowfrq[i,2]] <- lowfrq[i,3]
        }
      }

      if(ignoregappositions == TRUE){
        gapfound <- apply(dfseq, 1, function(x) paste(names(which(x=="-"))))
        gapposition <- unique(unlist(gapfound))
        for (i in gapposition){
          dfseq[[i]] <- "-"
        }
      }else if(ignoregappositions == FALSE | missing(ignoregappositions)){
        dfseq <- dfseq
      }else{
        dfseq <- dfseq
      }

      dfseq2 <- unite(dfseq, col='seq', c(names(dfseq[1:ncol(dfseq)])), sep='')
      readnames <- dss2df(seq)
      dfseq2$names <- readnames$names
      seq <- DNAStringSet(c(dfseq2$seq))
      seq@ranges@NAMES <- c(dfseq2$names)

      ##sub-sampling after
      if(samplingfirst == TRUE | missing(samsize)){
        seq <- seq
      }else if(samplingfirst == FALSE | missing(samplingfirst)){
        if(samsize < length(seq)){
          samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
          seq <- seq[names(seq) %in% samreads]
        }else{
          seq <- seq
        }
      }else{
        if(samsize < length(seq)){
          samreads <- sample(seq@ranges@NAMES, samsize, replace = FALSE)
          seq <- seq[names(seq) %in% samreads]
        }else{
          seq <- seq
        }
      }

      ##vqs analysis
      hap <- QSutils::Collapse(seq)
      if(length(hap$hseqs) <= 2){
        hapre <- hap
        names(hapre)[2] <- "seqs"
      }else{
        hapcor <- CorrectGapsAndNs(hap$hseqs[2:length(hap$hseqs)],hap$hseqs[[1]])
        hapcor <- c(hap$hseqs[1],hapcor)
        hapre <- QSutils::Recollapse(hapcor,hap$nr)
      }

      depth <- length(seq)
      nsingleton <- sum(hapre$nr == 1)
      pctsingleton <- nsingleton*100/depth

      storeout <- data.frame(label, pct, pctsingleton)
      store <- rbind(store, storeout)

      if(pctsingleton <= pctsing) {
        break
      }
    }
    if(missing(label)){
      return(pct)
    }else{
      return(store)
    }
  }
}
