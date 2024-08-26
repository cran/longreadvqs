#' Comparing operational taxonomic unit (OTU) by k-means clustering between samples
#'
#' @description
#' Pools noise-minimized down-sampled read samples and compares their diversity based on operational taxonomic unit (OTU) classified by k-means clustering of single nucleotide variant (SNV) distance. This function is a subset of "vqscompare" function.
#'
#' @param samplelist List of samples, i.e., name of resulting objects from "vqsassess" or "vqscustompct" functions, for example list(BC1, BC2, BC3).
#' @param kmeans.n Number of operational taxonomic units (OTUs) needed from k-means clustering on multidimensional scale (MDS) of all samples' pairwise genetic distance.
#'
#' @return Comparative table of OTU diversity metrics between listed samples calculated from consensus sequence of each OTU by QSutils package
#' @export
#'
#' @import dplyr
#' @importFrom stats as.dist
#' @importFrom stats cmdscale
#' @importFrom stats kmeans
#' @importFrom ape as.DNAbin
#' @importFrom ape dist.dna
#' @importFrom plyr rbind.fill
#'
#' @examples
#' ## Locate input FASTA files-----------------------------------------------------------------------
#' sample1filepath <- system.file("extdata", "s1.fasta", package = "longreadvqs")
#' sample2filepath <- system.file("extdata", "s2.fasta", package = "longreadvqs")
#'
#' ## Prepare data for viral quasispecies comparison between two samples-----------------------------
#' sample1 <- vqsassess(sample1filepath, pct = 10, samsize = 20, label = "sample1")
#' sample2 <- vqsassess(sample2filepath, pct = 10, samsize = 20, label = "sample2")
#'
#' ## Compare OTU (4 clusters) diversity metrics between two samples---------------------------------
#' otucompare(samplelist = list(sample1, sample2), kmeans.n = 4)
#'
#' @name otucompare

otucompare <- function(samplelist = list(BC1, BC2, BC3), kmeans.n = 20){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)), sapply(sapply(LETTERS, function(x) paste0(x, LETTERS)), function(x) paste0(x, LETTERS)))

  allhap <- do.call(c, purrr::map(samplelist, "hapre"))
  comb <- plyr::rbind.fill(purrr::map(samplelist, "snvhap"))
  bar <- comb[,c("group", "value", "label", "hap")]
  lab <- as.character(do.call(c, purrr::map(samplelist, "lab")))

  datall <- dss2df(allhap)
  sumall <- datall %>% group_by(seq) %>% summarize(names = paste0(unique(names), collapse = ',')) %>% arrange(desc(str_length(names)), names)
  sumall$newgroup <- LETTERS702[1:nrow(sumall)]
  sumgroup <- sumall[,2:3] %>% mutate(names = strsplit(as.character(names), ",")) %>% unnest(names)
  colnames(sumgroup) <- c("hap", "newgroup")

  sumallseq <- data.frame(str_split_fixed(sumall$seq, "", max(nchar(sumall$seq))))
  colnames(sumallseq) <- paste0(1:ncol(sumallseq))
  sumallsnv <- sumallseq[vapply(sumallseq, function(x) length(unique(x)) > 1, logical(1L))]
  sumallsnv <- unite(sumallsnv, col='seq', c(names(sumallsnv[1:ncol(sumallsnv)])), sep='')
  sumallsnv <- DNAStringSet(c(sumallsnv$seq))
  sumallsnv@ranges@NAMES <- c(sumall$newgroup)
  sumallsnv <- as.DNAbin(sumallsnv)
  distdnasnv <- dist.dna(sumallsnv, pairwise.deletion = FALSE, model = "raw")
  dist_sumallsnv <- as.dist(distdnasnv)
  mds_sumallsnv <- cmdscale(dist_sumallsnv) %>% as_tibble(.name_repair = 'unique')
  colnames(mds_sumallsnv) <- c("Dim.1", "Dim.2")
  clust <- kmeans(mds_sumallsnv, kmeans.n)$cluster %>% as.factor()
  mds_sumallsnv <- mds_sumallsnv %>% mutate(otu = clust)
  mds_sumallsnv$newgroup <- attr(dist_sumallsnv, "Labels")
  sumgroup <- inner_join(sumgroup, mds_sumallsnv, by = "newgroup")

  bar <- merge(x=bar,y=sumgroup,by="hap",all.x=TRUE)
  bar_hap <- bar %>% group_by(newgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))
  bar_hap$sortgrp <- LETTERS702[1:nrow(bar_hap)]
  bar_hap <- merge(x=bar,y=bar_hap,by="newgroup",all.x=TRUE) %>% arrange(desc(nhap))
  bar_hap <- bar_hap %>% group_by(sortgrp) %>% mutate(mean = mean(nhap)) %>%
    ungroup() %>% arrange(desc(mean)) %>% mutate(sortgrp = factor(sortgrp, levels = unique(sortgrp)))
  bar_hap <- bar_hap %>% group_by(label) %>% mutate(pct = value/sum(value)) %>% ungroup()
  bar_hap_ori <- bar_hap[, c("sortgrp", "hap", "value")] %>% dplyr::rename("hapgroup" = "sortgrp")

  bar_otu <- bar %>% group_by(otu) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))
  bar_otu$sortgrp <- LETTERS702[1:nrow(bar_otu)]
  bar_otu <- merge(x=bar,y=bar_otu,by="otu",all.x=TRUE) %>% arrange(desc(nhap))
  bar_otu_ori <- bar_otu[, c("sortgrp", "hap")] %>% dplyr::rename("otugroup" = "sortgrp")

  hap_otu <- inner_join(bar_hap_ori, bar_otu_ori, by = "hap")
  hap_otu <- inner_join(hap_otu, bar_hap[c('hap', 'sortgrp')], by = "hap") %>% dplyr::rename("id" = "hap", "hapgroup2" = "sortgrp")
  hap_otu[c('group', 'number')] <- str_split_fixed(hap_otu$id, '_', 2)

  fullseq_hap_otu <- inner_join(datall, hap_otu[, c("id", "group", "hapgroup", "otugroup", "value")], by = c("names" = "id"))
  uniq.otu <- unique(fullseq_hap_otu$otugroup) %>% as.list()
  d.otu <- data.frame()

  for (i in uniq.otu) {
    eaotu <- subset(fullseq_hap_otu, otugroup == i)
    mat.eaotu <- data.frame(str_split_fixed(eaotu$seq, "", max(nchar(eaotu$seq)))) %>% as.matrix()
    con.eaotu <- seqinr::consensus(mat.eaotu)
    d.otu.i <- data.frame(consensus = paste(con.eaotu, collapse = ''), otugroup = i)
    d.otu <- rbind(d.otu,d.otu.i)
  }

  sumotugrp <- fullseq_hap_otu %>% group_by(group, otugroup) %>% summarise(nreads = sum(value)) %>% as.data.frame()
  sumotugrp <- inner_join(sumotugrp, d.otu, by = "otugroup")
  sumotugrp <- sumotugrp %>% mutate(width = nchar(consensus))

  uniq.sample <- unique(sumotugrp$group) %>% as.list()
  d.otudiversity <- data.frame()

  for (i in uniq.sample) {
    easample <- subset(sumotugrp, group == i) %>% arrange(desc(nreads))
    hapresample <- list()
    hapresample$nr <- easample$nreads
    hapresample$seqs <- DNAStringSet(c(easample$consensus))
    hapresample$seqs@ranges@NAMES <- c(easample$otugroup)

    label <- i
    depth <- sum(hapresample$nr)
    OTUs <- length(hapresample$seqs)
    polymorph <- SegSites(hapresample$seqs)
    mutations <- TotalMutations(hapresample$seqs)
    shannon <- Shannon(hapresample$nr)
    norm_shannon <- NormShannon(hapresample$nr)
    gini_simpson <- QSutils::GiniSimpson(hapresample$nr)
    ##Functional diversity
    dst <- DNA.dist(hapresample$seqs,model="raw")
    #Incidence-based (count)
    FAD <- FAD(dst)
    Mfe <- MutationFreq(dst)
    Pie <- NucleotideDiversity(dst)
    #Abundance-based (frequency)
    nm <- nmismatch(pairwiseAlignment(hapresample$seqs,hapresample$seqs[1]))
    Mfm <- MutationFreqVar(nm,hapresample$nr,len=width(hapresample$seqs)[1])
    Pim <- NucleotideDiversity(dst,hapresample$nr)
    d.otudiversity.i <- data.frame(label, depth, OTUs, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    d.otudiversity <- rbind(d.otudiversity,d.otudiversity.i)
  }

  res <- d.otudiversity
  return(res)
}
