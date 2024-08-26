#' Comparing viral quasispecies profile and operational taxonomic unit (OTU) classified by k-means clustering between samples
#'
#' @description
#' Pools noise-minimized down-sampled read samples and compares their diversity by 1) viral quasispecies profile (haplotype and metrics from QSutils package), 2) operational taxonomic unit (OTU) classified by k-means clustering of single nucleotide variant (SNV) distance, and 3) visualization of different comparative method, i.e., haplotype, OTU, phylogenetic tree, MDS plot. Such comparisons can also be performed at the amino acid level (protein haplotype and single amino acid variation (SAV) group).
#'
#' @param samplelist List of samples, i.e., name of resulting objects from "vqsassess" or "vqscustompct" functions, for example list(BC1, BC2, BC3).
#' @param lab_name Name of variable or type of sample for instance "barcode", "sample", "dpi", or "isolate" (optional).
#' @param kmeans.n Number of operational taxonomic units (OTUs) and single amino acid variation (SAV) groups needed from k-means clustering on multidimensional scale (MDS) of all samples' pairwise SNV and SAV distances.
#' @param showhap.n Number of largest haplotypes (default = 30) labeled in the top five OTUs' MDS plot (optional).
#' @param proteincoding Translate gene or protein-coding reads into amino acid sequences and regroup them into protein haplotypes and single amino acid variation (SAV) groups, which are comparable to haplotypes and OTUs at the nucleotide level, respectively (optional). If not specified or if proteincoding = FALSE, gene translation and downstream analyses will not be performed (default).
#' @param removestopcodon Remove the last amino acid (expected to be a stop codon) from translated amino acid sequences before further analysis (optional). If not specified or if removestopcodon = FALSE, the last amino acid will not be removed (default).
#'
#' @return List of 1) "hapdiv": comparative table of viral quasispecies diversity metrics between listed samples calculated by QSutils package, 2) "otudiv": comparative table of operational taxonomic unit (OTU) diversity metrics between listed samples calculated from consensus sequence of each OTU (similar to "otucompare" function's output), 3) "sumsnv_hap_otu": frequency and SNV profile (by position in the alignment) of all haplotypes and OTUs, 4) "fullseq": complete read sequence of all haplotypes, 5) "fulldata": complete read sequence of all haplotypes in every sample with frequency and OTU classification, 6) "summaryplot": visualization of viral quasispecies comparison between samples including  6.1) "happlot": proportion of haplotypes (top left), 6.2) "otuplot": proportion of OTUs (bottom left), 6.3) multidimensional scale (MDS) plots (right) of k-means OTU ("top5otumds": 5 largest groups with major haplotypes labeled and "allotumds": all groups), 7) "aadiv": comparative table of viral quasispecies diversity metrics between listed samples based on translated reads calculated by QSutils package (similar to one of "AAcompare" function's outputs), 8) "savgrpdiv": comparative table of single amino acid (SAV) group diversity metrics between listed samples calculated from consensus amino acid sequence of each SAV group (similar to one of "AAcompare" function's outputs), 9) "sumsav_phap_savgrp": frequency and SAV profile (by position in the alignment) of all protein haplotypes and SAV groups, 10) "fullseq_aa": complete amino acid sequence of all protein haplotypes, 11) "fulldata_aa": complete amino acid sequence of all protein haplotypes in every sample with frequency and SAV group classification, 12) "summaryplot_aa": visualization of viral quasispecies comparison between samples based on translated reads including  12.1) "phapplot": proportion of protein haplotypes (top left), 12.2) "savgrpplot": proportion of SAV groups (bottom left), 12.3) multidimensional scale (MDS) plots (right) of k-means SAV group ("top5savgrpmds": 5 largest SAV groups with major protein haplotypes labeled and "allsavgrpmds": all SAV groups), 7) to 12) will be generated only when proteincoding = TRUE
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
#' @importFrom ape as.AAbin
#' @importFrom ape dist.dna
#' @importFrom ape dist.aa
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

utils::globalVariables(c('element_text','guides','hapgroup2','newgroup','nhap','nreads','hap','otu','otugroup','AAgroup','savgroup','hapgroup','nphap','pct','position_fill','sortgrp','xlab','geom_text','scale_color_manual'))

vqscompare <- function(samplelist = list(BC1, BC2, BC3), lab_name = "sample name", kmeans.n = 20, showhap.n = 30, proteincoding = FALSE, removestopcodon = FALSE){
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

  if(proteincoding == TRUE){
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
  }

  sumallseq <- data.frame(str_split_fixed(sumall$seq, "", max(nchar(sumall$seq))))
  colnames(sumallseq) <- paste0(1:ncol(sumallseq))
  sumallsnv <- sumallseq[vapply(sumallseq, function(x) length(unique(x)) > 1, logical(1L))]
  combsnv <- sumallsnv
  combsnv$hap <- sumall$names
  combsnv <- combsnv %>% separate_longer_delim(hap, delim = ",")
  sumallsnv <- unite(sumallsnv, col='seq', c(names(sumallsnv[1:ncol(sumallsnv)])), sep='')
  sumallsnv <- DNAStringSet(c(sumallsnv$seq))
  sumallsnv@ranges@NAMES <- c(sumall$newgroup)
  sumallsnv <- as.DNAbin(sumallsnv)
  distdnasnv <- dist.dna(sumallsnv, pairwise.deletion = FALSE, model = "raw")
  dist_sumallsnv <- as.dist(distdnasnv)
  mds_sumallsnv <- suppressMessages({cmdscale(dist_sumallsnv) %>% as_tibble(.name_repair = 'unique')})
  colnames(mds_sumallsnv) <- c("Dim.1", "Dim.2")
  clust <- kmeans(mds_sumallsnv, kmeans.n)$cluster %>% as.factor()
  mds_sumallsnv <- mds_sumallsnv %>% mutate(otu = clust)
  mds_sumallsnv$newgroup <- attr(dist_sumallsnv, "Labels")
  sumgroup <- inner_join(sumgroup, mds_sumallsnv, by = "newgroup")

  if(proteincoding == TRUE){
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
    sumgroup <- inner_join(sumgroup, sumgroup_aa, by = "hap")
  }

  bar <- merge(x=bar,y=sumgroup,by="hap",all.x=TRUE)
  bar_hap <-suppressMessages({bar %>% group_by(newgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
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

  if(proteincoding == TRUE){
    bar_aa <- suppressMessages({bar %>% group_by(AAgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
    bar_aa$sortgrp <- letters702[1:nrow(bar_aa)]
    bar_aa <- merge(x=bar,y=bar_aa,by="AAgroup",all.x=TRUE) %>% arrange(desc(nhap))
    bar_aa_ori <- bar_aa[, c("sortgrp", "hap")] %>% dplyr::rename("AAgroup" = "sortgrp")
    bar_aa_sum <- aggregate(value ~ sortgrp + label, data = bar_aa, FUN = sum)
    bar_aa_sum <- bar_aa_sum %>% group_by(sortgrp) %>% mutate(mean = mean(value)) %>%
      ungroup() %>% arrange(desc(mean)) %>% mutate(sortgrp = factor(sortgrp, levels = unique(sortgrp)))
    bar_aa_sum <- bar_aa_sum %>% group_by(label) %>% mutate(pct = value/sum(value)) %>% ungroup()
    bar_aa_sum$sortgrp[bar_aa_sum$value == 1] <- NA
    bar_aa_sum$sortgrp <- as.character(bar_aa_sum$sortgrp)
    bar_aa_sum$sortgrp <- factor(bar_aa_sum$sortgrp, levels=unique(bar_aa_sum$sortgrp)[order(nchar(unique(bar_aa_sum$sortgrp)), unique(bar_aa_sum$sortgrp))])
    colaa <- sample(color, length(unique(bar_aa_sum$sortgrp)), replace = TRUE)
    aaprop <- ggplot(bar_aa_sum, aes(fill = sortgrp, x=factor(label, levels = lab), y = value)) + geom_bar(colour = "black", position = "fill", stat = "identity") +
      scale_fill_manual(values = colaa) + ggtitle("Protein haplotypes") + theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("proportion") + xlab(lab_name) +
      geom_text(aes(label = ifelse(pct > .05, as.character(sortgrp), "")), position = position_fill(vjust = .5))
  }

  bar_otu <- suppressMessages({bar %>% group_by(otu) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
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

  if(proteincoding == TRUE){
    bar_savgroup <- suppressMessages({bar %>% group_by(savgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
    bar_savgroup$sortgrp <- LETTERS702[1:nrow(bar_savgroup)]
    bar_savgroup <- merge(x=bar,y=bar_savgroup,by="savgroup",all.x=TRUE) %>% arrange(desc(nhap))
    bar_savgroup_ori <- bar_savgroup[, c("sortgrp", "hap")] %>% dplyr::rename("savgroupsort" = "sortgrp")
    bar_savgroup_sum <- aggregate(value ~ sortgrp + label, data = bar_savgroup, FUN = sum)
    bar_savgroup_sum <- bar_savgroup_sum %>% group_by(sortgrp) %>% mutate(mean = mean(value)) %>%
      ungroup() %>% arrange(desc(mean)) %>% mutate(sortgrp = factor(sortgrp, levels = unique(sortgrp)))
    bar_savgroup_sum <- bar_savgroup_sum %>% group_by(label) %>% mutate(pct = value/sum(value)) %>% ungroup()
    bar_savgroup_sum$sortgrp <- as.character(bar_savgroup_sum$sortgrp)
    bar_savgroup_sum$sortgrp <- factor(bar_savgroup_sum$sortgrp, levels=unique(bar_savgroup_sum$sortgrp)[order(nchar(unique(bar_savgroup_sum$sortgrp)), unique(bar_savgroup_sum$sortgrp))])
    colsavgroup <- sample(color, length(unique(bar_savgroup_sum$sortgrp)), replace = TRUE)
    savgroupprop <- ggplot(bar_savgroup_sum, aes(fill = sortgrp, x=factor(label, levels = lab), y = value)) + geom_bar(colour = "black", position = "fill", stat = "identity") +
      scale_fill_manual(values = colsavgroup) + ggtitle("SAV groups") + theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("proportion") + xlab(lab_name) +
      geom_text(aes(label = ifelse(pct > .03, as.character(sortgrp), "")), position = position_fill(vjust = .5))
  }

  if(proteincoding == TRUE){
    comb <- combsnv %>% select(hap, everything())
    allsnv <- subset(bar_hap, select = -c(otu, AAgroup, savgroup)) %>% left_join(comb, by="hap") %>% dplyr::rename("hapgroup" = "sortgrp")
    otulist <- bar_otu[, c("hap", "sortgrp")] %>% dplyr::rename("otugroup" = "sortgrp")
    allsnv <- left_join(allsnv, otulist, by = "hap")
    aalist <- bar_aa[, c("hap", "sortgrp")] %>% dplyr::rename("AAgroup" = "sortgrp")
    allsnv <- left_join(allsnv, aalist, by = "hap")
    savlist <- bar_savgroup[, c("hap", "sortgrp")] %>% dplyr::rename("savgroup" = "sortgrp")
    allsnv <- left_join(allsnv, savlist, by = "hap")
    combsav <- combsav %>% select(hap, everything())
    allsav <- left_join(allsnv[, c("hap", "value", "hapgroup", "otugroup", "AAgroup", "savgroup")], combsav, by = "hap")
  } else {
    comb <- combsnv %>% select(hap, everything())
    allsnv <- subset(bar_hap, select = -c(otu)) %>% left_join(comb, by="hap") %>% dplyr::rename("hapgroup" = "sortgrp")
    otulist <- bar_otu[, c("hap", "sortgrp")] %>% dplyr::rename("otugroup" = "sortgrp")
    allsnv <- left_join(allsnv, otulist, by = "hap")
  }

  fullnt_otu <- allsnv[, c("hap", "nhap", "otugroup")] %>% left_join(datall[, c("seq", "names")], by=c('hap'='names'))
  fullnt2_otu <- left_join(fullnt_otu, bar_hap_ori[, c("hap", "hapgroup")], by ="hap")
  fullnt2_otu <- fullnt2_otu[, c("otugroup", "hapgroup", "nhap", "seq")] %>% distinct() %>% arrange(str_length(otugroup), otugroup, desc(nhap))
  fullntseq_otu <- data.frame(str_split_fixed(fullnt2_otu$seq, "", max(nchar(fullnt2_otu$seq))))
  colnames(fullntseq_otu) <- paste0(1:ncol(fullntseq_otu))
  fullntsnvs_otu <- fullntseq_otu[vapply(fullntseq_otu, function(x) length(unique(x)) > 1, logical(1L))]
  compsnv_otu <- cbind(fullnt2_otu[, c("otugroup", "hapgroup", "nhap")], fullntsnvs_otu)

  if(proteincoding == TRUE){
    compsavcol <- suppressMessages({subset(allsav, select = -c(hap, hapgroup, otugroup, savgroup)) %>% group_by(AAgroup) %>% summarise(nphap = sum(value)) %>% arrange(desc(nphap))})
    compsavgroup <- distinct(subset(allsav, select = -c(hap, value, hapgroup, otugroup)))
    compsavgroup <- left_join(compsavcol, compsavgroup, by = "AAgroup") %>% arrange(str_length(savgroup), savgroup, desc(nphap))
    compsavgroup <- compsavgroup %>% select(savgroup, AAgroup, nphap, everything())
  }

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
  allplot <- plot_grid(propcomb, pmdscomb, ncol = 2, rel_widths = c(1,1))

  if(proteincoding == TRUE){
    mdsdat_aa <- inner_join(allsav[, c("hap", "savgroup", "AAgroup")], bar_savgroup[, c("hap", "Dim.1_aa", "Dim.2_aa")], by = "hap")
    mdsdat_aa <- unique(subset(mdsdat_aa, select=-c(hap)))
    pmds_aa <- suppressWarnings({ggpubr::ggscatter(mdsdat_aa, x = "Dim.1_aa", y = "Dim.2_aa", color = "savgroup", palette = colsavgroup, size =1, ellipse = TRUE, repel = TRUE) +
        ggtitle("All SAV groups") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))})
    mdsdat2_aa <- mdsdat_aa

    if(missing(showhap.n)){
      mdsdat2_aa$AAgroup[!(mdsdat2_aa$AAgroup %in% letters702[1:30])] <- ""
    }else{
      mdsdat2_aa$AAgroup[!(mdsdat2_aa$AAgroup %in% letters702[1:showhap.n])] <- ""
    }

    pmds2_aa <- suppressWarnings({ggpubr::ggscatter(subset(mdsdat2_aa, savgroup %in% c("A", "B", "C", "D", "E")), x = "Dim.1_aa", y = "Dim.2_aa", color = "savgroup", palette = colsavgroup, size =1, ellipse = TRUE, repel = TRUE, label = "AAgroup") +
        ggtitle("Top5 SAV groups") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))})

    propcomb_aa <- plot_grid(aaprop, savgroupprop, nrow = 2)
    pmdscomb_aa <- suppressWarnings({plot_grid(pmds2_aa, pmds_aa, nrow = 2)})
    allplot_aa <- plot_grid(propcomb_aa, pmdscomb_aa, ncol = 2, rel_widths = c(1,1))
  }

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

  sumotugrp <- suppressMessages({fullseq_hap_otu %>% group_by(group, otugroup) %>% summarise(nreads = sum(value)) %>% as.data.frame()})
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

  if(proteincoding == TRUE){
    aa_sav <- inner_join(bar_aa_ori, bar_savgroup_ori, by = "hap")
    aa_sav <- inner_join(aa_sav, bar_aa[c('hap', 'value')], by = "hap") %>% dplyr::rename("id" = "hap", "savgroup" = "savgroupsort")
    aa_sav[c('group', 'number')] <- str_split_fixed(aa_sav$id, '_', 2)
    fullseq_aa_sav <- inner_join(datall_aa, aa_sav[, c("id", "group", "AAgroup", "savgroup", "value")], by = c("names" = "id"))
    fullaa2_sav <- suppressMessages({fullseq_aa_sav %>% group_by(AAgroup) %>% summarise(nhap = sum(value)) %>% arrange(desc(nhap))})
    tempfullaa2 <- fullseq_aa_sav[c('savgroup', 'AAgroup', 'seq')] %>% distinct()
    fullaa2_sav <- inner_join(fullaa2_sav, tempfullaa2, by = "AAgroup") %>% arrange(str_length(savgroup), savgroup, desc(nhap))
    fullaa2_sav <- fullaa2_sav[c('savgroup', 'AAgroup', 'nhap', 'seq')]
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
  }

  if(proteincoding == TRUE){
    list <- list("summaryplot" = allplot, "happlot" = happrop, "otuplot" = otuprop, "top5otumds" = pmds2, "allotumds" = pmds,
                 "summaryplot_aa" = allplot_aa, "phapplot" = aaprop, "savgrpplot" = savgroupprop, "top5savgrpmds" = pmds2_aa, "allsavgrpmds" = pmds_aa,
                 "hapdiv" = d.hapdiversity, "otudiv" = d.otudiversity, "aadiv" = d.AAdiversity, "savgrpdiv" = d.savdiversity,
                 "sumsnv_hap_otu" = compsnv_otu, "sumsav_phap_savgrp" = compsavgroup,
                 "fullseq" = data.frame(fullnt2_otu[, c("hapgroup", "seq")]), "fullseq_aa" = data.frame(fullaa2_sav[, c("AAgroup", "seq")]),
                 "fulldata" = fullseq_hap_otu[,-1], "fulldata_aa" = fullseq_aa_sav[,-1])
  } else {
    list <- list("summaryplot" = allplot, "happlot" = happrop, "otuplot" = otuprop, "top5otumds" = pmds2, "allotumds" = pmds,
                 "hapdiv" = d.hapdiversity, "otudiv" = d.otudiversity,
                 "sumsnv_hap_otu" = compsnv_otu, "fullseq" = data.frame(fullnt2_otu[, c("hapgroup", "seq")]), "fulldata" = fullseq_hap_otu[,-1])
  }

  return(list)
}
