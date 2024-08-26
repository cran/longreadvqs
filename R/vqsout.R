#' Exporting viral quasispecies profile comparison results
#'
#' @description
#' Writes out resulting objects from "vqscompare" function as tables (TSV files) and alignment (FASTA file) to the working directory.
#'
#' @param vqscompare.obj A resulting object from "vqscompare" function.
#' @param directory Path to desired directory (location) for output files. If it is not specified, the directory will be the current working directory.
#'
#' @return TSV files of viral quasispecies profile comparison results and FASTA file of unique haplotype alignment.
#' @export
#'
#' @importFrom utils write.table
#'
#' @examples
#' ## Locate input FASTA files-----------------------------------------------------------------------
#' sample1filepath <- system.file("extdata", "s1.fasta", package = "longreadvqs")
#' sample2filepath <- system.file("extdata", "s2.fasta", package = "longreadvqs")
#'
#' ## Prepare data for viral quasispecies comparison between two samples-----------------------------
#' set.seed(123)
#' sample1 <- vqsassess(sample1filepath, pct = 0, samsize = 50, label = "sample1")
#' sample2 <- vqsassess(sample2filepath, pct = 0, samsize = 50, label = "sample2")
#'
#' ## Compare viral quasispecies and OTU (4 clusters) diversity between two samples------------------
#' comp <- vqscompare(samplelist = list(sample1, sample2),
#'            lab_name = "Sample", kmeans.n = 4, showhap.n = 5)
#'
#' ## Export Key outputs from "vqscompare" function--------------------------------------------------
#' notrun <- vqsout(comp, directory = tempdir())
#'
#' @name vqsout

vqsout <- function(vqscompare.obj, directory = "path/to/directory"){
  if(missing(directory)) {
    write.table(vqscompare.obj$hapdiv, file = "hap_div.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$otudiv, file = "otu_div.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$sumsnv_hap_otu, file = "snv_hap_otu.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$fulldata, file = "fulldata.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
    writeFasta<-function(data, filename){
      fastaLines = c()
      for (rowNum in 1:nrow(data)){
        fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"hapgroup"], sep = "")))
        fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
      }
      fileConn<-file(filename)
      writeLines(fastaLines, fileConn)
      close(fileConn)
    }
    writeFasta(vqscompare.obj$fullseq, "fullseq.fasta")
    if(length(vqscompare.obj) > 10){
      write.table(vqscompare.obj$aadiv, file = "phap_div.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$savgrpdiv, file = "savgrp_div.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$sumsav_phap_savgrp, file = "sav_phap_savgrp.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$fulldata_aa, file = "fulldata_aa.tsv", sep = '\t', row.names = FALSE, quote = FALSE)
      writeFasta<-function(data, filename){
        fastaLines = c()
        for (rowNum in 1:nrow(data)){
          fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"AAgroup"], sep = "")))
          fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
        }
        fileConn<-file(filename)
        writeLines(fastaLines, fileConn)
        close(fileConn)
      }
      writeFasta(vqscompare.obj$fullseq_aa, "fullseq_aa.fasta")
    }
  }else{
    write.table(vqscompare.obj$hapdiv, file = file.path(directory, paste("hap_div.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$otudiv, file = file.path(directory, paste("otu_div.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$sumsnv_hap_otu, file = file.path(directory, paste("snv_hap_otu.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(vqscompare.obj$fulldata, file = file.path(directory, paste("fulldata.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
    writeFasta<-function(data, filename){
      fastaLines = c()
      for (rowNum in 1:nrow(data)){
        fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"hapgroup"], sep = "")))
        fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
      }
      fileConn<-file(filename)
      writeLines(fastaLines, fileConn)
      close(fileConn)
    }
    writeFasta(vqscompare.obj$fullseq, file.path(directory, paste("fullseq.fasta")))
    if(length(vqscompare.obj) > 10){
      write.table(vqscompare.obj$aadiv, file = file.path(directory, paste("phap_div.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$savgrpdiv, file = file.path(directory, paste("savgrp_div.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$sumsav_phap_savgrp, file = file.path(directory, paste("sav_phap_savgrp.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
      write.table(vqscompare.obj$fulldata_aa, file = file.path(directory, paste("fulldata_aa.tsv")), sep = '\t', row.names = FALSE, quote = FALSE)
      writeFasta<-function(data, filename){
        fastaLines = c()
        for (rowNum in 1:nrow(data)){
          fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"AAgroup"], sep = "")))
          fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
        }
        fileConn<-file(filename)
        writeLines(fastaLines, fileConn)
        close(fileConn)
      }
      writeFasta(vqscompare.obj$fullseq_aa, file.path(directory, paste("fullseq_aa.fasta")))
    }
  }
}
