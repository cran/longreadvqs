#' Plotting single nucleotide variant (SNV) frequency in read alignment across different samples
#'
#' @description
#' Compares single nucleotide variant (SNV) profile between error-minimized down-sampled read samples using cowplot's "plot_grid" function. The resulting plot may help evaluating the optimal cut-off percentage of low frequency nucleotide base used in "vqsassess", "vqscustompct", or "vqssub" functions.
#'
#' @param samplelist List of samples, i.e., name of resulting objects from "vqsassess" or "vqscustompct" functions, for example list(BC1, BC2, BC3).
#' @param ncol Number of columns for multiple plots (see cowplot's "plot_grid" function)
#'
#' @return Comparative plot of SNV frequency in read alignment across different samples
#' @export
#'
#' @import scales
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
#'
#' @examples
#' ## Locate input FASTA files-----------------------------------------------------------------------
#' sample1filepath <- system.file("extdata", "s1.fasta", package = "longreadvqs")
#' sample2filepath <- system.file("extdata", "s2.fasta", package = "longreadvqs")
#'
#' ## Prepare data for viral quasispecies comparison between two samples-----------------------------
#' sample1 <- vqsassess(sample1filepath, pct = 10, label = "sample1")
#' sample2 <- vqsassess(sample2filepath, pct = 10, label = "sample2")
#'
#' ## Compare SNV profile between two listed samples-------------------------------------------------
#' snvcompare(samplelist = list(sample1, sample2), ncol = 1)
#'
#' @name snvcompare

utils::globalVariables(c('BC1','BC2','BC3','geom_col','plot_grid','theme'))

snvcompare <- function(samplelist = list(BC1, BC2, BC3), ncol = 1, barwidth = 1){
  samplelist = samplelist
  list.plot <- purrr::map(samplelist, "snv") %>% as.list()
  list.plot2 <- list()
  for(i in 1:length(list.plot)) {
    list.plot2[[i]] <- list.plot[[i]] + theme(legend.position = "none")
  }
  legend <- suppressWarnings({cowplot::get_legend(list.plot[[1]])})
  snvplot0 <- suppressWarnings({cowplot::plot_grid(plotlist = list.plot2, ncol = ncol, align = "v")})
  snvplot <- suppressWarnings({cowplot::plot_grid(snvplot0, legend, ncol = 2, rel_widths = c(1, .1))})
  return(snvplot)
}
