#' Sequencing error minimization with customized % cut-off at particular nucleotide region, read down-sampling, and data preparation for viral quasispecies comparison
#'
#' @description
#' Minimizes potential long-read sequencing error based on the specified cut-off percentages of low frequency nucleotide base and down-samples read for further comparison with other samples. In this function, the cut-off percentage can be specifically adjusted for different ranges of nucleotide positions which is very useful when sequencing error heavily occurs in a particular part of reads. The output of this function is a list of several objects representing diversity of each sample that must be used as an input for other functions such as "snvcompare" or "vqscompare".
#'
#' @param fasta Input as a read alignment in FASTA format
#' @param method Sequencing error minimization methods that replace low frequency nucleotide base (less than the "pct" cut-off) with consensus base of that position ("conbase": default) or with base of the dominant haplotype ("domhapbase").
#' @param samplingfirst Downsampling before (TRUE) or after (FALSE: default) the error minimization.
#' @param pct Percent cut-off defining low frequency nucleotide base that will be replaced (must be specified).
#' @param brkpos Ranges of nucleotide positions with different % cut-off specified in "lspct" for example c("1:50","51:1112") meaning that the first and the second ranges are nucleotide positions 1 to 50 and 51 to 1112, respectively.
#' @param lspct List of customized % cut-off applied to nucleotide ranges set in "brkpos" for example c(15,8) meaning that 15% and 8% cut-offs will be applied to the first and the second ranges, respectively.
#' @param gappct The percent cut-off particularly specified for gap (-). If it is not specified or less than "pct", "gappct" will be equal to "pct" (default).
#' @param ignoregappositions Replace all nucleotides in the positions in the alignment containing gap(s) with gap. This will make such positions no longer single nucleotide variant (SNV). The default is "FALSE".
#' @param samsize Sample size (number of reads) after down-sampling. If it is not specified or more than number of reads in the original alignment, down-sampling will not be performed (default).
#' @param label String within quotation marks indicating name of read alignment (optional). Please don't use underscore (_) in the label.
#'
#' @return list of 1) "dat": viral quasispecies diversity metrics calculated by QSutils package (similar to "vqssub" function's output), 2) "snvhap": SNV profile of each haplotype with frequency and new label for "vqscompare" function, 3) "snv": plot of SNV frequency for "snvcompare" function, 4) "hapre": DNAStringSet of read alignment of each haplotype for "vqscompare" function, 5) "lab": name of sample or read alignment
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
#' @importFrom reshape2 melt
#' @importFrom magrittr set_colnames
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import QSutils
#'
#' @examples
#' ## Locate input FASTA file------------------------------------------------------------------------
#' fastafilepath <- system.file("extdata", "badend.fasta", package = "longreadvqs")
#'
#' ## Prepare data for viral quasispecies comparison using 10% cut-off across all positions----------
#' nocustom <- vqsassess(fastafilepath, pct = 10, label = "nocustom")
#'
#' ## Prepare data using 10% cut-off for the first 74 positions and 30% cut-off for the rest---------
#' custom <- vqscustompct(fastafilepath, pct = 10,
#'                        brkpos = c("1:74","75:84"), lspct = c(10,30), label = "custom")
#'
#' ## Use "snvcompare" function to check whether SNV profile looks better or not---------------------
#' snvcompare(samplelist = list(nocustom, custom), ncol = 1)
#'
#' @name vqscustompct

utils::globalVariables(c('newcol','aes','base','geom_bar','ggplot','ggtitle','group','percent','position','scale_fill_manual','value','ylab','set_colnames'))

vqscustompct <- function(fasta, method= c("conbase", "domhapbase"), samplingfirst = FALSE, pct = 8, brkpos = c("1:50","51:width(seq[1])"), lspct = c(15,8), gappct = 50, ignoregappositions = TRUE, samsize = 100, label = "sample"){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)), sapply(sapply(LETTERS, function(x) paste0(x, LETTERS)), function(x) paste0(x, LETTERS)))
  cols <- c("-"  = "#ffffbf", "A" = "#d7191c", "G" ="#fdae61", "C" = "#abdda4", "T" = "#2b83ba")
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

  if(grepl("_", label) == FALSE){
    label <- label
  }else if(grepl("_", label) == TRUE){
    stop("Please rename the label without underscore symbol (_)")
  }else{
    label <- label
  }

  method <- match.arg(method)

  if(method == "conbase" | missing(method)){
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

    if(missing(brkpos)){
      stop("Please specify a list of ranges of nucleotide positions (brkpos)")
    }else if(missing(lspct)){
      stop("Please specify a list of different % cut-off for low frequency SNV (lspct)")
    }else if(length(brkpos) != length(lspct)){
      stop("Please specify brkpos and lspct as lists with the same length")
    }else{
      datalist <- list()
      for (i in 1:length(lspct)){
        pos <- eval(parse(text = brkpos[i]))
        datalist[[i]] <- data.frame(apply(frqpc[pos,], 1, function(x) paste(names(which(x[1:5] < lspct[i] & x[1:5] > 0)), collapse = ",")))
      }
      lowfrq <- do.call(rbind, datalist) %>% set_colnames("lowfrq")
    }

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

    seqcor <- CorrectGapsAndNs(seq[1:length(seq)],hap$hseqs[[1]])
    dfhap <- dss2df(seqcor)
    dfhap <- data.frame(str_split_fixed(dfhap$seq, "", max(nchar(dfhap$seq))))
    mathap <- as.matrix(dfhap)
    frqhap <- seqinr::consensus(mathap, method = "profile")
    frqhap <- as.data.frame(frqhap)
    frqhappc <- 100 * sweep(frqhap, 2, colSums(frqhap), `/`)
    colnames(frqhappc) <- c(1:ncol(frqhappc))
    frqhap2 <- melt(as.matrix(frqhappc), value.name = "percent", varnames=c('base', 'position'))
    frqhap2$percent[frqhap2$percent == 100] <- 0

    bar <- data.frame(group=LETTERS702[1:length(hapre$nr)], value=hapre$nr)
    dfhapre <- dss2df(hapre$seqs)
    dfhapre <- data.frame(str_split_fixed(dfhapre$seq, "", max(nchar(dfhapre$seq))))
    colnames(dfhapre) <- paste0(1:ncol(dfhapre))
    dfhapre <- dfhapre[,colnames(dfhapre) %in% unique(subset(frqhap2, percent != 0)$position), drop=FALSE]
    dfhapre$group <- LETTERS702[1:length(hapre$nr)]

    depth <- length(seq)
    haplotypes <- length(hapre$seqs)
    polymorph <- SegSites(hapre$seqs)
    mutations <- TotalMutations(hapre$seqs)

    shannon <- Shannon(hapre$nr)
    norm_shannon <- NormShannon(hapre$nr)
    gini_simpson <- QSutils::GiniSimpson(hapre$nr)

    ##Functional diversity
    dst <- DNA.dist(hapre$seqs,model="raw")
    #Incidence-based (count)
    FAD <- FAD(dst)
    Mfe <- MutationFreq(dst)
    Pie <- NucleotideDiversity(dst)
    #Abundance-based (frequency)
    nm <- nmismatch(pairwiseAlignment(hapre$seqs,hapre$seqs[1]))
    Mfm <- MutationFreqVar(nm,hapre$nr,len=width(hapre$seqs)[1])
    Pim <- NucleotideDiversity(dst,hapre$nr)

    nsingleton <- sum(hapre$nr == 1)
    pctsingleton <- nsingleton*100/depth

    dat <- data.frame(label, method, samplingwhen, pct, fulldepth, depth, haplotypes, nsingleton, pctsingleton, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    snv <- ggplot(frqhap2, aes(fill = base, x = position, y = percent)) + geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(values = cols) + ylab("proportion") + ggtitle(label)
    snvhap <- inner_join(bar,dfhapre, by = "group")
    snvhap_newlab <- snvhap
    snvhap_newlab$label <- label
    snvhap_newlab$hap <- 1:nrow(snvhap_newlab)
    snvhap_newlab$hap <- paste(label, snvhap_newlab$hap, sep='_')
    snvhap_newlab <- snvhap_newlab %>% select(label, hap, group, value, everything())
    bar <- snvhap_newlab[,c("group", "value", "label", "hap")]
    hapre_newlab <- hapre
    hapre_newlab$seqs@ranges@NAMES <- bar$hap
    list <- list("dat" = dat, "snv" = snv, "hapre" = hapre_newlab$seqs, "snvhap"= snvhap_newlab, lab = label)
    return(list)
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

    if(missing(brkpos)){
      stop("Please specify a list of ranges of nucleotide positions (brkpos)")
    }else if(missing(lspct)){
      stop("Please specify a list of different % cut-off for low frequency SNV (lspct)")
    }else if(length(brkpos) != length(lspct)){
      stop("Please specify brkpos and lspct as lists with the same length")
    }else{
      datalist <- list()
      for (i in 1:length(lspct)){
        pos <- eval(parse(text = brkpos[i]))
        datalist[[i]] <- data.frame(apply(frqpc[pos,], 1, function(x) paste(names(which(x[1:5] < lspct[i] & x[1:5] > 0)), collapse = ",")))
      }
      lowfrq <- do.call(rbind, datalist) %>% set_colnames("lowfrq")
    }

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
      message(paste0("Dominant haplotype is ", format(round(domhap$nr[1]*100/sum(domhap$nr), 2), nsmall = 2), "% of total reads"))
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

    seqcor <- CorrectGapsAndNs(seq[1:length(seq)],hap$hseqs[[1]])
    dfhap <- dss2df(seqcor)
    dfhap <- data.frame(str_split_fixed(dfhap$seq, "", max(nchar(dfhap$seq))))
    mathap <- as.matrix(dfhap)
    frqhap <- seqinr::consensus(mathap, method = "profile")
    frqhap <- as.data.frame(frqhap)
    frqhappc <- 100 * sweep(frqhap, 2, colSums(frqhap), `/`)
    colnames(frqhappc) <- c(1:ncol(frqhappc))
    frqhap2 <- melt(as.matrix(frqhappc), value.name = "percent", varnames=c('base', 'position'))
    frqhap2$percent[frqhap2$percent == 100] <- 0

    bar <- data.frame(group=LETTERS702[1:length(hapre$nr)], value=hapre$nr)
    dfhapre <- dss2df(hapre$seqs)
    dfhapre <- data.frame(str_split_fixed(dfhapre$seq, "", max(nchar(dfhapre$seq))))
    colnames(dfhapre) <- paste0(1:ncol(dfhapre))
    dfhapre <- dfhapre[,colnames(dfhapre) %in% unique(subset(frqhap2, percent != 0)$position), drop=FALSE]
    dfhapre$group <- LETTERS702[1:length(hapre$nr)]

    depth <- length(seq)
    haplotypes <- length(hapre$seqs)
    polymorph <- SegSites(hapre$seqs)
    mutations <- TotalMutations(hapre$seqs)

    shannon <- Shannon(hapre$nr)
    norm_shannon <- NormShannon(hapre$nr)
    gini_simpson <- QSutils::GiniSimpson(hapre$nr)

    ##Functional diversity
    dst <- DNA.dist(hapre$seqs,model="raw")
    #Incidence-based (count)
    FAD <- FAD(dst)
    Mfe <- MutationFreq(dst)
    Pie <- NucleotideDiversity(dst)
    #Abundance-based (frequency)
    nm <- nmismatch(pairwiseAlignment(hapre$seqs,hapre$seqs[1]))
    Mfm <- MutationFreqVar(nm,hapre$nr,len=width(hapre$seqs)[1])
    Pim <- NucleotideDiversity(dst,hapre$nr)

    nsingleton <- sum(hapre$nr == 1)
    pctsingleton <- nsingleton*100/depth

    dat <- data.frame(label, method, samplingwhen, pct, fulldepth, depth, haplotypes, nsingleton, pctsingleton, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    snv <- ggplot(frqhap2, aes(fill = base, x = position, y = percent)) + geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(values = cols) + ylab("proportion") + ggtitle(label)
    snvhap <- inner_join(bar,dfhapre, by = "group")
    snvhap_newlab <- snvhap
    snvhap_newlab$label <- label
    snvhap_newlab$hap <- 1:nrow(snvhap_newlab)
    snvhap_newlab$hap <- paste(label, snvhap_newlab$hap, sep='_')
    snvhap_newlab <- snvhap_newlab %>% select(label, hap, group, value, everything())
    bar <- snvhap_newlab[,c("group", "value", "label", "hap")]
    hapre_newlab <- hapre
    hapre_newlab$seqs@ranges@NAMES <- bar$hap
    list <- list("dat" = dat, "snv" = snv, "hapre" = hapre_newlab$seqs, "snvhap"= snvhap_newlab, lab = label)
    return(list)
  }
}
