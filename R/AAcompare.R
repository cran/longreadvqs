#' Comparing viral quasispecies diversity metrics at amino acid level
#'
#' @description
#' Pools noise-minimized down-sampled read samples and compares their diversity metrices based on protein haplotype and single amino acid variation (SAV) group that is classified by k-means clustering of SAV distance. This function is a subset of "vqscompare" function.
#'
#' @param samplelist List of samples, i.e., name of resulting objects from "vqsassess" or "vqscustompct" functions, for example list(BC1, BC2, BC3).
#' @param kmeans.n Number of single amino acid variation (SAV) groups needed from k-means clustering on multidimensional scale (MDS) of all samples' pairwise SAV distance.
#' @param removestopcodon Remove the last amino acid (expected to be a stop codon) from translated amino acid sequences before further analysis (optional). If not specified or if removestopcodon = FALSE, the last amino acid will not be removed (default).
#'
#' @return List of 1) "aadiv": comparative table of viral quasispecies diversity metrics between listed samples based on translated reads calculated by QSutils package, and 2) "savgrpdiv": comparative table of single amino acid (SAV) group diversity metrics between listed samples calculated from consensus amino acid sequence of each SAV group
#' @export
#'
#' @import dplyr
#' @importFrom stats as.dist
#' @importFrom stats cmdscale
#' @importFrom stats kmeans
#' @importFrom ape as.AAbin
#' @importFrom ape dist.aa
#' @importFrom plyr rbind.fill
#'
#' @examples
#' ## Locate input FASTA files-----------------------------------------------------------------------
#' sample1filepath <- system.file("extdata", "s1.fasta", package = "longreadvqs")
#' sample2filepath <- system.file("extdata", "s2.fasta", package = "longreadvqs")
#'
#' ## Prepare data for viral quasispecies comparison between two samples-----------------------------
#' set.seed(123)
#' sample1 <- vqsassess(sample1filepath, pct = 5, samsize = 50, label = "sample1")
#' sample2 <- vqsassess(sample2filepath, pct = 5, samsize = 50, label = "sample2")
#'
#' ## Compare protein haplotype and SAV group (4 clusters) diversity metrics between two samples-----
#' AAcompare(samplelist = list(sample1, sample2), kmeans.n = 4)
#'
#' @name AAcompare

AAcompare <- function(samplelist = list(BC1, BC2, BC3), kmeans.n = 20, removestopcodon = FALSE){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)), sapply(sapply(LETTERS, function(x) paste0(x, LETTERS)), function(x) paste0(x, LETTERS)))

  allhap <- do.call(c, purrr::map(samplelist, "hapre"))
  comb <- plyr::rbind.fill(purrr::map(samplelist, "snvhap"))
  bar <- comb[,c("group", "value", "label", "hap")]
  lab <- as.character(do.call(c, purrr::map(samplelist, "lab")))

  letters702 <- c(letters, sapply(letters, function(x) paste0(x, letters)), sapply(sapply(letters, function(x) paste0(x, letters)), function(x) paste0(x, letters)))
  getCodons <- function(myAln) {
    seqs <- as.character(myAln)
    len <- width(myAln)[1]
    starts <- seq(from=1, to=len, by=3)
    ends <- starts + 2
    myViews <- lapply(myAln, function(x) {
      Views(x, starts, ends)
    })
    myCodons <- lapply(myViews, function(x) {
      as.character(DNAStringSet(x))
    })
    myCodons
  }
  translateCodons <- function(myCodons, unknownCodonTranslatesTo="-") {
    ## make new genetic code
    gapCodon <- "-"
    names(gapCodon) <- "---"
    my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)

    ## translate the codons
    pep <- my_GENETIC_CODE[myCodons]

    ## check for codons that were not possible to translate, e.g. frameshift codons
    if (sum(is.na(pep))>0) {
      pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
    }

    ## prep for output
    pep <- paste(pep, collapse="")
    return(pep)
  }
  translateGappedAln <- function(myAln, unknownCodonTranslatesTo="-") {
    myCodons <- getCodons(myAln)
    myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
    return(myAAaln)
  }
  allhap_aa <- translateGappedAln(allhap, unknownCodonTranslatesTo="X")
  datall_aa <- dss2df(allhap_aa)
  if(removestopcodon == TRUE | missing(removestopcodon)){
    datall_aa$seq = substr(datall_aa$seq,1,nchar(datall_aa$seq)-1)
    datall_aa$width = datall_aa$width - 1
  }
  sumall_aa <- datall_aa %>% group_by(seq) %>% summarize(names = paste0(unique(names), collapse = ',')) %>% arrange(desc(str_length(names)), names)
  sumall_aa$newgroup <- LETTERS702[1:nrow(sumall_aa)]
  sumgroup_aa <- sumall_aa[,2:3] %>% mutate(names = strsplit(as.character(names), ",")) %>% unnest(names)
  colnames(sumgroup_aa) <- c("hap", "AAgroup")

  sumallseq_aa <- data.frame(str_split_fixed(sumall_aa$seq, "", max(nchar(sumall_aa$seq))))
  colnames(sumallseq_aa) <- paste0(1:ncol(sumallseq_aa))
  sumallsav <- sumallseq_aa[vapply(sumallseq_aa, function(x) length(unique(x)) > 1, logical(1L))]
  combsav <- sumallsav
  combsav$hap <- sumall_aa$names
  combsav <- combsav %>% separate_longer_delim(hap, delim = ",")
  sumallsav <- unite(sumallsav, col='seq', c(names(sumallsav[1:ncol(sumallsav)])), sep='')
  sumallsav <- AAStringSet(c(sumallsav$seq))
  sumallsav@ranges@NAMES <- c(sumall_aa$newgroup)
  savlength <- width(sumallsav)[1]
  sumallsav <- as.AAbin(sumallsav)
  sumallsav <- as.matrix(sumallsav)
  distaasav <- ape::dist.aa(sumallsav)
  dist_sumallsav <- as.dist(distaasav)
  mds_sumallsav <- suppressMessages({cmdscale(dist_sumallsav) %>% as_tibble(.name_repair = 'unique')})
  colnames(mds_sumallsav) <- c("Dim.1_aa", "Dim.2_aa")
  clust_aa <- kmeans(mds_sumallsav, kmeans.n)$cluster %>% as.factor()
  mds_sumallsav <- mds_sumallsav %>% mutate(savgroup = clust_aa)
  mds_sumallsav$AAgroup <- attr(dist_sumallsav, "Labels")
  sumgroup_aa <- inner_join(sumgroup_aa, mds_sumallsav, by = "AAgroup")

  bar <- merge(x=bar,y=sumgroup_aa,by="hap",all.x=TRUE)
  bar_aa <- suppressMessages({bar %>% group_by(AAgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
  bar_aa$sortgrp <- letters702[1:nrow(bar_aa)]
  bar_aa <- merge(x=bar,y=bar_aa,by="AAgroup",all.x=TRUE) %>% arrange(desc(nhap))
  bar_aa_ori <- bar_aa[, c("sortgrp", "hap")] %>% dplyr::rename("AAgroup" = "sortgrp")

  bar_savgroup <- suppressMessages({bar %>% group_by(savgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
  bar_savgroup$sortgrp <- LETTERS702[1:nrow(bar_savgroup)]
  bar_savgroup <- merge(x=bar,y=bar_savgroup,by="savgroup",all.x=TRUE) %>% arrange(desc(nhap))
  bar_savgroup_ori <- bar_savgroup[, c("sortgrp", "hap")] %>% dplyr::rename("savgroupsort" = "sortgrp")

  aa_sav <- inner_join(bar_aa_ori, bar_savgroup_ori, by = "hap")
  aa_sav <- inner_join(aa_sav, bar_aa[c('hap', 'value')], by = "hap") %>% dplyr::rename("id" = "hap", "savgroup" = "savgroupsort")
  aa_sav[c('group', 'number')] <- str_split_fixed(aa_sav$id, '_', 2)
  fullseq_aa_sav <- inner_join(datall_aa, aa_sav[, c("id", "group", "AAgroup", "savgroup", "value")], by = c("names" = "id"))
  uniq.sav <- unique(fullseq_aa_sav$savgroup) %>% as.list()
  d.sav <- data.frame()

  for (i in uniq.sav) {
    easav <- subset(fullseq_aa_sav, savgroup == i)
    mat.easav <- data.frame(str_split_fixed(easav$seq, "", max(nchar(easav$seq)))) %>% as.matrix()
    con.easav <- seqinr::consensus(mat.easav)
    d.sav.i <- data.frame(consensus = paste(con.easav, collapse = ''), savgroup = i)
    d.sav <- rbind(d.sav,d.sav.i)
  }

  sumsavgrp <- suppressMessages({fullseq_aa_sav %>% group_by(group, savgroup) %>% summarise(nreads = sum(value)) %>% as.data.frame()})
  sumsavgrp <- inner_join(sumsavgrp, d.sav, by = "savgroup")
  sumsavgrp <- sumsavgrp %>% mutate(width = nchar(consensus))

  uniq.sample_sav <- unique(sumsavgrp$group) %>% as.list()
  d.savdiversity <- data.frame()

  for (i in uniq.sample_sav) {
    easample <- subset(sumsavgrp, group == i) %>% arrange(desc(nreads))
    hapresample <- list()
    hapresample$nr <- easample$nreads
    hapresample$seqs <- AAStringSet(c(easample$consensus))
    hapresample$seqs@ranges@NAMES <- c(easample$savgroup)

    label <- i
    depth <- sum(hapresample$nr)
    SAVgroups <- length(hapresample$seqs)
    polymorph <- SegSites(hapresample$seqs)
    mutations <- TotalMutations(hapresample$seqs)
    shannon <- Shannon(hapresample$nr)
    norm_shannon <- NormShannon(hapresample$nr)
    gini_simpson <- QSutils::GiniSimpson(hapresample$nr)
    ##Functional diversity
    AA.dist <- function(seqs){
      if(!is(seqs, "AAStringSet"))
        stop("The input object must be AAStringSet \n")
      # Convert the alignment into AAbin object
      strm <- as.AAbin(seqs)
      strm <- as.matrix(strm)
      # Compute the matrix of distances
      dst <- ape::dist.aa(strm)
      # Convert NAs into 0
      dst[is.na(dst)] <- 0
      return(dst)
    }
    dst <- AA.dist(hapresample$seqs)
    #Incidence-based (count)
    FAD <- FAD(dst)
    Mfe <- MutationFreq(dst)
    Pie <- NucleotideDiversity(dst)
    #Abundance-based (frequency)
    nm <- nmismatch(pairwiseAlignment(hapresample$seqs,hapresample$seqs[1]))
    Mfm <- MutationFreqVar(nm,hapresample$nr,len=width(hapresample$seqs)[1])
    Pim <- NucleotideDiversity(dst,hapresample$nr)
    nsingleton <- sum(hapresample$nr == 1)
    pctsingleton <- nsingleton*100/depth
    d.savdiversity.i <- data.frame(label, depth, SAVgroups, nsingleton, pctsingleton, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    d.savdiversity <- rbind(d.savdiversity,d.savdiversity.i)
  }

  sumAAgrp <- suppressMessages({fullseq_aa_sav %>% group_by(group, AAgroup) %>% summarise(nreads = sum(value)) %>% as.data.frame()})
  fstemp <- fullseq_aa_sav[c('AAgroup', 'seq')] %>% distinct()
  sumAAgrp <- inner_join(sumAAgrp, fstemp, by = "AAgroup")
  sumAAgrp <- sumAAgrp %>% mutate(width = nchar(seq))

  uniq.sample_aa <- unique(sumAAgrp$group) %>% as.list()
  d.AAdiversity <- data.frame()

  for (i in uniq.sample_aa) {
    easample <- subset(sumAAgrp, group == i) %>% arrange(desc(nreads))
    hapresample <- list()
    hapresample$nr <- easample$nreads
    hapresample$seqs <- AAStringSet(c(easample$seq))
    hapresample$seqs@ranges@NAMES <- c(easample$AAgroup)

    label <- i
    depth <- sum(hapresample$nr)
    protein_haplotypes <- length(hapresample$seqs)
    polymorph <- SegSites(hapresample$seqs)
    mutations <- TotalMutations(hapresample$seqs)
    shannon <- Shannon(hapresample$nr)
    norm_shannon <- NormShannon(hapresample$nr)
    gini_simpson <- QSutils::GiniSimpson(hapresample$nr)
    ##Functional diversity
    AA.dist <- function(seqs){
      if(!is(seqs, "AAStringSet"))
        stop("The input object must be AAStringSet \n")
      # Convert the alignment into AAbin object
      strm <- as.AAbin(seqs)
      strm <- as.matrix(strm)
      # Compute the matrix of distances
      dst <- ape::dist.aa(strm)
      # Convert NAs into 0
      dst[is.na(dst)] <- 0
      return(dst)
    }
    dst <- AA.dist(hapresample$seqs)
    #Incidence-based (count)
    FAD <- FAD(dst)
    Mfe <- MutationFreq(dst)
    Pie <- NucleotideDiversity(dst)
    #Abundance-based (frequency)
    nm <- nmismatch(pairwiseAlignment(hapresample$seqs,hapresample$seqs[1]))
    Mfm <- MutationFreqVar(nm,hapresample$nr,len=width(hapresample$seqs)[1])
    Pim <- NucleotideDiversity(dst,hapresample$nr)
    nsingleton <- sum(hapresample$nr == 1)
    pctsingleton <- nsingleton*100/depth
    d.AAdiversity.i <- data.frame(label, depth, protein_haplotypes, nsingleton, pctsingleton, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    d.AAdiversity <- rbind(d.AAdiversity,d.AAdiversity.i)
  }

  list <- list("aadiv" = d.AAdiversity, "savgrpdiv" = d.savdiversity)
  return(list)
}
