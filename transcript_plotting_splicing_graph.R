## script to make transcript blocks above splicing graph

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Gviz")

## script to make transcript blocks above splicing graph
library(SplicingGraphs)
library(Gviz)
library(AnnotationDbi)
library(GenomicFeatures)
library(igraph)
library(tidyverse)

txdb = loadDb(file = 'txdb.gencode34.sqlite')
isActiveSeq(txdb)[1:25]
isActiveSeq(txdb)[-match("chr1", names(isActiveSeq(txdb)))] <- FALSE #just make chr1 active
sg <- SplicingGraphs(txdb)
edges_by_gene <- sgedgesByGene(sg)

#obtain splice graph info from gene of interest
geneID <- "ENSG00000228826.3"
edg_gene <- as.data.frame(sgedges(sg[geneID]))

#obtain how many transcripts in gene
gene_edges <- as.data.frame(edges_by_gene[[geneID]])
exon_edges <- gene_edges %>% filter(ex_or_in %in% "ex")
tx_list <- c()
for(i in 1:length(gene_edges$tx_id)){
  tx_list <- c(tx_list, gene_edges$tx_id[[i]])
}
tx_list <- unique(tx_list)

#splicing graph plot
par(mar=c(7,15,7,1), bty="n")
plot(0,0, type="n",xaxt="n", yaxt="n", xlab=geneID, ylab="",
     xlim=c(1,length(sgnodes(sg[geneID]))), ylim=c(-(nrow(edg_gene)/10),length(tx_list)+2))

x <- 1:length(sgnodes(sg[geneID]))
y <- rep(1, length(sgnodes(sg[geneID])))

points(x, y+0*2, cex=3, col="orange")
text(x,y,labels=c("R",1:(length(sgnodes(sg[geneID]))-2),"L"))

iArrows <- igraph:::igraph.Arrows

for(i in 1:nrow(edg_gene)) {
  if(edg_gene$ex_or_in[i]=="ex"){
    iArrows(as.numeric(edg_gene$from[i])+1, 1, as.numeric(edg_gene$to[i])+1, 1,
            h.lwd=0.25, sh.lwd=1, sh.col="darkgreen",
            width=2, size=0.1)
  }
  if(edg_gene$ex_or_in[i]=="in"){
    iArrows(as.numeric(edg_gene$from[i])+1, 1, as.numeric(edg_gene$to[i])+1, 1,
            h.lwd=0.25, sh.lwd=1, sh.col="orange",
            width=2, size=0.1, curve=0.3 - (i %% 2), sh.lty=2)
  }
  if(edg_gene$ex_or_in[i]=="" & edg_gene$from[i]=="R"){
    iArrows(1, 1, as.numeric(edg_gene$to[i])+1, 1,
            h.lwd=0.25, sh.lwd=1, sh.col="orange",
            width=2, size=0.1, curve=0.3 - (i %% 2), sh.lty=2,)
  }
  if(edg_gene$ex_or_in[i]=="" & edg_gene$to[i]=="L"){
    iArrows(as.numeric(edg_gene$from[i])+1, 1, length(sgnodes(sg[geneID])), 1,
            h.lwd=0.25, sh.lwd=1, sh.col="orange",
            width=2, size=0.1, curve=0.3 - (i %% 2), sh.lty=2)
  }
}


#transcript plot
count = 0
for(i in 1:length(tx_list)){
  x0 <- c()
  x1 <- c()
  y <- rep(2+count, length(grep(tx_list[i], exon_edges$tx_id)))
  for(j in grep(tx_list[i], exon_edges$tx_id)){
    x0 <- c(x0, as.numeric(exon_edges$from[j])+1)
    x1 <- c(x1, as.numeric(exon_edges$to[j])+1)
  }
  text(x=-1, y=y, tx_list[i], cex=0.8, xpd = TRUE)
  segments(x0 = x0, y0 = y, x1 = x1, y1 = y, col="darkgreen", lwd=2)
  count = count + 1
}
