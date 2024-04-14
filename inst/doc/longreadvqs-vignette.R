## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(longreadvqs)
library(ggplot2)

#Load sample 1.
sample1 <- system.file("extdata", "sample1.fasta", package = "longreadvqs")

## ----fig0, warning=FALSE, fig.height = 3, fig.width = 7.2, fig.align = "center"----
#Check which % cut-off could effectively minimize errors by assessing % singleton haplotypes.
x <- pctopt(sample1, pctsing = 0, label = "sample1")
ggplot(x, aes(x=pct, y=pctsingleton)) + geom_line() + geom_point()

## ----cmd0, warning=FALSE, message=FALSE---------------------------------------
#VQS diversity metrics of error minimized (10% cut-off) read alignment
vqssub(sample1, pct = 10, label = "sample1")

## ----cmd1, warning=FALSE, message=FALSE---------------------------------------
#Load samples 2 and 3.
sample2 <- system.file("extdata", "sample2.fasta", package = "longreadvqs")
sample3 <- system.file("extdata", "sample3.fasta", package = "longreadvqs")

#Error minimization (10% cut-off) and down-sampling (depth of 300)
a <- vqssub(sample1, pct = 10, samsize = 300, label = "sample1")
b <- vqssub(sample2, pct = 10, samsize = 300, label = "sample2")
c <- vqssub(sample3, pct = 10, samsize = 300, label = "sample3")

#Compare VQS diversity across three samples after error minimization and down-sampling.
rbind(a, b, c)

## ----cmd2, warning=FALSE, message=FALSE---------------------------------------
#Randomly down-sample from depth of 655 to 300, and 100 (10 iterations for each sample size).
set.seed(123)
c_full <- vqssub(sample3, pct = 10, label = "full_depth") #non-sampled alignment
c_300 <- vqsresub(sample3, iter = 10, pct = 10, 
                  samsize = 300, label = "depth_300") #down-sample to depth of 300
c_100 <- vqsresub(sample3, iter = 10, pct = 10,
                  samsize = 100, label = "depth_100") #down-sample to depth of 100

all_c <- rbind(c_full, c_300, c_100) #combine data

#Load packages for visualization.
library(cowplot)

all_c$depth <- as.character(all_c$depth)
depthorder <- c("655", "300", "100")

#Plot box-plots of three VQS metrics variation through different sampling depths.
haplotypes <- ggplot(all_c, aes(x = factor(depth, level=depthorder), y = haplotypes, 
	group = interaction(depth, label), fill = label)) + geom_boxplot() + 
  xlab("depth") + theme(legend.position = "none")
shannon <- ggplot(all_c, aes(x = factor(depth, level=depthorder), y = shannon, 
	group = interaction(depth, label), fill = label)) + geom_boxplot() + 
  xlab("depth") + theme(legend.position = "none")
gini_simpson <- ggplot(all_c, aes(x = factor(depth, level=depthorder), y = gini_simpson, 
	group = interaction(depth, label), fill = label)) + geom_boxplot() + 
  xlab("depth") + theme(legend.position = "none")

## ----fig1, fig.height = 3, fig.width = 7.2, fig.align = "center"--------------
plot_grid(haplotypes, shannon, gini_simpson, nrow = 1)

## ----cmd3, warning=FALSE, message=FALSE---------------------------------------
#Generate complete VQS data for further comparison for each sample.
s1 <- vqsassess(sample1, pct = 10, samsize = 300, label = "sample1")
s2 <- vqsassess(sample2, pct = 10, samsize = 300, label = "sample2")
s3 <- vqsassess(sample3, pct = 10, samsize = 300, label = "sample3")

## ----fig2, fig.height = 6, fig.width = 7.2, fig.align = "center"--------------
#Compare SNV profile between three samples.
snvcompare(samplelist = list(s1, s2, s3), ncol = 1)

## ----cmd4, warning=FALSE, message=FALSE---------------------------------------
#Load sample 4.
sample4 <- system.file("extdata", "mock.fasta", package = "longreadvqs")
s4 <- vqsassess(sample4, pct = 10, samsize = 300, label = "sample4")

## ----fig3, warning=FALSE, fig.height = 2, fig.width = 7.2, fig.align = "center"----
#SNV profile of sample 4
s4$snv

## ----cmd5, warning=FALSE, message=FALSE---------------------------------------
#Use "vqscustompct" function to increase % cut-off after position 972 to 30%.
s4_fix <- vqscustompct(sample4, pct = 10, brkpos = c("1:971","972:982"), lspct = c(10,30), 
                       samsize = 300, label = "sample4")

## ----fig4, warning=FALSE, fig.height = 2, fig.width = 7.2, fig.align = "center"----
#SNV profile of sample 4 after % cut-off adjustment
s4_fix$snv

## ----cmd6, warning=FALSE, message=FALSE---------------------------------------
#List outputs from "vqsassess" or "vqscustompct" that we want to compare into "vqscompare" function.
#Set the number of new OTU groups based on k-means clustering to 10 groups (kmeans.n = 10).
set.seed(1234)
comp <- vqscompare(samplelist = list(s1, s2, s3, s4_fix),
                   lab_name = "Sample", kmeans.n = 10)

## ----fig5, warning=FALSE, fig.height = 6, fig.width = 7.2, fig.align = "center"----
#The most important output of the "vqscompare" function is the summary plot.
comp$summaryplot

