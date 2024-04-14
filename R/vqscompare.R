#' Comparing viral quasispecies profile and operational taxonomic unit (OTU) classified by k-means clustering between samples
#'
#' @description
#' Pools error-minimized down-sampled read samples and compares their diversity by 1) viral quasispecies profile (haplotype and metrics from QSutils package), 2) operational taxonomic unit (OTU) classified by k-means clustering of single nucleotide variant (SNV) distance, and 3) visualization of different comparative method, i.e., haplotype, OTU, phylogenetic tree, MDS plot.
#'
#' @param samplelist List of samples, i.e., name of resulting objects from "vqsassess" or "vqscustompct" functions, for example list(BC1, BC2, BC3).
#' @param lab_name Name of variable or type of sample for instance "barcode", "sample", "dpi", or "isolate" (optional).
#' @param kmeans.n Number of clusters or operational taxonomic units (OTUs) needed from k-means clustering on multidimensional scale (MDS) of all samples' pairwise SNV distance.
#' @param showhap.n Number of largest haplotypes (default = 30) labeled in the top five OTUs' MDS plot (optional).
#'
#' @return list of 1) "hapdiv": comparative table of viral quasispecies diversity metrics between listed samples calculated by QSutils package, 2) "otudiv": comparative table of OTU diversity metrics between listed samples calculated from consensus sequence of each OTU (similar to "otucompare" function's output), 3) "sumsnv_hap": frequency and SNV profile (by position in the alignment) of haplotypes that are not singleton (number of reads > 1), 4) "sumsnv_otu": frequency and SNV profile of all haplotypes grouped into different operational taxonomic unit (OTU), 5) "fullseq": complete read sequence of haplotypes that are not singleton, 6) "fulldata": complete read sequence of all haplotypes in every sample with frequency and OTU classification, 7) "summaryplot": visualization of viral quasispecies comparison between samples including 7.1) "happlot": proportion of haplotypes (top left), 7.2) "otuplot": proportion of OTUs (bottom left), and 7.3) multidimensional scale (MDS) plots (right) of k-means OTU ("top5otumds": 5 largest groups with major haplotypes labeled and "allotumds": all groups)
#' @export
#'
#' @import RColorBrewer
#' @import dplyr
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom stats as.dist
#' @importFrom stats cmdscale
#' @importFrom stats kmeans
#' @importFrom ape as.DNAbin
#' @importFrom ape bionjs
#' @importFrom ape dist.dna
#' @importFrom ggpubr ggscatter
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
#' ## Compare viral quasispecies and OTU (4 clusters) diversity between two samples------------------
#' out <- vqscompare(samplelist = list(sample1, sample2),
#'            lab_name = "Sample", kmeans.n = 4, showhap.n = 5)
#' out$summaryplot
#'
#' @name vqscompare

utils::globalVariables(c('element_text','guides','hapgroup2','newgroup','nhap','nreads','hap','otu','otugroup','pct','position_fill','sortgrp','xlab','geom_text','scale_color_manual'))

vqscompare <- function(samplelist = list(BC1, BC2, BC3), lab_name = "sample name", kmeans.n = 20, showhap.n = 30){
  dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)), sapply(sapply(LETTERS, function(x) paste0(x, LETTERS)), function(x) paste0(x, LETTERS)))

  allhap <- do.call(c, purrr::map(samplelist, "hapre"))
  comb <- plyr::rbind.fill(purrr::map(samplelist, "snvhap"))
  bar <- comb[,c("group", "value", "label", "hap")]
  lab <- as.character(do.call(c, purrr::map(samplelist, "lab")))
  d.hapdiversity <- plyr::rbind.fill(purrr::map(samplelist, "dat"))

  if(missing(lab_name)){
    lab_name = "sample name"
    }

  datall <- dss2df(allhap)
  sumall <- datall %>% group_by(seq) %>% summarize(names = paste0(unique(names), collapse = ',')) %>% arrange(desc(str_length(names)), names)
  sumall$newgroup <- LETTERS702[1:nrow(sumall)]
  sumgroup <- sumall[,2:3] %>% mutate(names = strsplit(as.character(names), ",")) %>% unnest(names)
  colnames(sumgroup) <- c("hap", "newgroup")
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

  sumallseq <- data.frame(str_split_fixed(sumall$seq, "", max(nchar(sumall$seq))))
  colnames(sumallseq) <- paste0(1:ncol(sumallseq))
  sumallsnv <- sumallseq[vapply(sumallseq, function(x) length(unique(x)) > 1, logical(1L))]
  sumallsnv <- unite(sumallsnv, col='seq', c(names(sumallsnv[1:ncol(sumallsnv)])), sep='')
  sumallsnv <- DNAStringSet(c(sumallsnv$seq))
  sumallsnv@ranges@NAMES <- c(sumall$newgroup)
  sumallsnv <- as.DNAbin(sumallsnv)
  distdnasnv <- dist.dna(sumallsnv, pairwise.deletion = FALSE, model = "raw")
  dist_sumallsnv <- as.dist(distdnasnv)
  mds_sumallsnv <- suppressWarnings({cmdscale(dist_sumallsnv) %>% as_tibble(.name_repair = 'unique')})
  colnames(mds_sumallsnv) <- c("Dim.1", "Dim.2")
  clust <- kmeans(mds_sumallsnv, kmeans.n)$cluster %>% as.factor()
  mds_sumallsnv <- mds_sumallsnv %>% mutate(otu = clust)
  mds_sumallsnv$newgroup <- attr(dist_sumallsnv, "Labels")
  sumgroup <- inner_join(sumgroup, mds_sumallsnv, by = "newgroup")

  bar <- merge(x=bar,y=sumgroup,by="hap",all.x=TRUE)
  bar_hap <- bar %>% group_by(newgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))
  bar_hap$sortgrp <- 1:nrow(bar_hap)
  bar_hap <- merge(x=bar,y=bar_hap,by="newgroup",all.x=TRUE) %>% arrange(desc(nhap))
  bar_hap <- bar_hap %>% group_by(sortgrp) %>% mutate(mean = mean(nhap)) %>%
    ungroup() %>% arrange(desc(mean)) %>% mutate(sortgrp = factor(sortgrp, levels = unique(sortgrp)))
  bar_hap <- bar_hap %>% group_by(label) %>% mutate(pct = value/sum(value)) %>% ungroup()
  bar_hap_ori <- bar_hap[, c("sortgrp", "hap", "value")] %>% dplyr::rename("hapgroup" = "sortgrp")
  bar_hap$sortgrp[bar_hap$value == 1] <- NA
  bar_hap$sortgrp <- as.character(bar_hap$sortgrp)
  bar_hap$sortgrp <- factor(bar_hap$sortgrp, levels=unique(bar_hap$sortgrp)[order(nchar(unique(bar_hap$sortgrp)), unique(bar_hap$sortgrp))])
  colhap <- sample(color, length(unique(bar_hap$sortgrp)), replace = TRUE)
  happrop <- ggplot(bar_hap, aes(fill = sortgrp, x=factor(label, levels = lab), y = value)) + geom_bar(colour = "black", position = "fill", stat = "identity") +
    scale_fill_manual(values = colhap) + ggtitle("Haplotypes") + theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("proportion") + xlab(lab_name) +
    geom_text(aes(label = ifelse(pct > .03, as.character(sortgrp), "")), position = position_fill(vjust = .5))

  bar_otu <- bar %>% group_by(otu) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))
  bar_otu$sortgrp <- LETTERS702[1:nrow(bar_otu)]
  bar_otu <- merge(x=bar,y=bar_otu,by="otu",all.x=TRUE) %>% arrange(desc(nhap))
  bar_otu_ori <- bar_otu[, c("sortgrp", "hap")] %>% dplyr::rename("otugroup" = "sortgrp")
  bar_otu_sum <- aggregate(value ~ sortgrp + label, data = bar_otu, FUN = sum)
  bar_otu_sum <- bar_otu_sum %>% group_by(sortgrp) %>% mutate(mean = mean(value)) %>%
    ungroup() %>% arrange(desc(mean)) %>% mutate(sortgrp = factor(sortgrp, levels = unique(sortgrp)))
  bar_otu_sum <- bar_otu_sum %>% group_by(label) %>% mutate(pct = value/sum(value)) %>% ungroup()
  bar_otu_sum$sortgrp <- as.character(bar_otu_sum$sortgrp)
  bar_otu_sum$sortgrp <- factor(bar_otu_sum$sortgrp, levels=unique(bar_otu_sum$sortgrp)[order(nchar(unique(bar_otu_sum$sortgrp)), unique(bar_otu_sum$sortgrp))])
  colotu <- sample(color, length(unique(bar_otu_sum$sortgrp)), replace = TRUE)
  otuprop <- ggplot(bar_otu_sum, aes(fill = sortgrp, x=factor(label, levels = lab), y = value)) + geom_bar(colour = "black", position = "fill", stat = "identity") +
    scale_fill_manual(values = colotu) + ggtitle("OTUs") + theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("proportion") + xlab(lab_name) +
    geom_text(aes(label = ifelse(pct > .03, as.character(sortgrp), "")), position = position_fill(vjust = .5))

  colorder <- as.character(sort(as.numeric(colnames(comb[,!(colnames(comb) %in% c("label", "hap", "group", "value"))]))))
  comb <- comb[, c("hap", colorder)]
  allsnv <- bar_hap %>% left_join(comb, by="hap") %>% dplyr::rename("hapgroup" = "sortgrp")
  otulist <- bar_otu[, c("hap", "sortgrp")] %>% dplyr::rename("otugroup" = "sortgrp")
  allsnv <- inner_join(allsnv, otulist, by = "hap")

  fullnt <- allsnv[, c("hap", "nhap", "hapgroup")] %>% left_join(datall[, c("seq", "names")], by=c('hap'='names'))
  fullnt2 <- fullnt[, c("hapgroup", "nhap", "seq")] %>% distinct() %>% drop_na()
  fullntseq <- data.frame(str_split_fixed(fullnt2$seq, "", max(nchar(fullnt2$seq))))
  colnames(fullntseq) <- paste0(1:ncol(fullntseq))
  fullntsnvs <- fullntseq[vapply(fullntseq, function(x) length(unique(x)) > 1, logical(1L))]
  compsnv <- cbind(fullnt2[, c("hapgroup", "nhap")], fullntsnvs)

  fullnt_otu <- allsnv[, c("hap", "nhap", "otugroup")] %>% left_join(datall[, c("seq", "names")], by=c('hap'='names'))
  fullnt2_otu <- fullnt_otu[, c("otugroup", "nhap", "seq")] %>% distinct() %>% arrange(str_length(otugroup), otugroup, desc(nhap))
  fullntseq_otu <- data.frame(str_split_fixed(fullnt2_otu$seq, "", max(nchar(fullnt2_otu$seq))))
  colnames(fullntseq_otu) <- paste0(1:ncol(fullntseq_otu))
  fullntsnvs_otu <- fullntseq_otu[vapply(fullntseq_otu, function(x) length(unique(x)) > 1, logical(1L))]
  compsnv_otu <- cbind(fullnt2_otu[, c("otugroup", "nhap")], fullntsnvs_otu)

  hap_otu <- inner_join(bar_hap_ori, bar_otu_ori, by = "hap")
  hap_otu <- inner_join(hap_otu, bar_hap[c('hap', 'sortgrp')], by = "hap") %>% dplyr::rename("id" = "hap", "hapgroup2" = "sortgrp")
  hap_otu[c('group', 'number')] <- str_split_fixed(hap_otu$id, '_', 2)

  mdsdat <- inner_join(allsnv[, c("hap", "otugroup", "hapgroup")], bar_otu[, c("hap", "Dim.1", "Dim.2")], by = "hap")
  mdsdat <- unique(subset(mdsdat, select=-c(hap)))
  pmds <- suppressWarnings({ggpubr::ggscatter(mdsdat, x = "Dim.1", y = "Dim.2", color = "otugroup", palette = colotu, size =1, ellipse = TRUE, repel = TRUE) +
      ggtitle("All OTUs") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))})
  mdsdat2 <- mdsdat
  mdsdat2$hapgroup <- as.numeric(mdsdat2$hapgroup)

  if(missing(showhap.n)){
    mdsdat2$hapgroup[mdsdat2$hapgroup > 30] <- ""
  }else{
    mdsdat2$hapgroup[mdsdat2$hapgroup > showhap.n] <- ""
  }

  pmds2 <- suppressWarnings({ggpubr::ggscatter(subset(mdsdat2, otugroup %in% c("A", "B", "C", "D", "E")), x = "Dim.1", y = "Dim.2", color = "otugroup", palette = colotu, size =1, ellipse = TRUE, repel = TRUE, label = "hapgroup") +
      ggtitle("Top5 OTUs") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))})

  propcomb <- plot_grid(happrop, otuprop, nrow = 2)
  pmdscomb <- suppressWarnings({plot_grid(pmds2, pmds, nrow = 2)})
  #clustcomb <- plot_grid(p3, pmdscomb, nrow = 2, rel_heights = c(2,1))
  allplot <- plot_grid(propcomb, pmdscomb, ncol = 2, rel_widths = c(1,1))

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
    nsingleton <- sum(hapresample$nr == 1)
    pctsingleton <- nsingleton*100/depth
    d.otudiversity.i <- data.frame(label, depth, OTUs, nsingleton, pctsingleton, polymorph, mutations, shannon, norm_shannon, gini_simpson, FAD, Mfe, Pie, Mfm, Pim)
    d.otudiversity <- rbind(d.otudiversity,d.otudiversity.i)
  }

  list <- list("summaryplot" = allplot, "happlot" = happrop, "otuplot" = otuprop, "top5otumds" = pmds2, "allotumds" = pmds, "hapdiv" = d.hapdiversity, "otudiv" = d.otudiversity, "sumsnv_hap" = compsnv, "sumsnv_otu" = compsnv_otu, "fullseq" = data.frame(fullnt2[, c("hapgroup", "seq")]), "fulldata" = fullseq_hap_otu[,-1])
  return(list)
}
