
library(rtracklayer)
library(GenomicFeatures)
library(SplicingGraphs)
library(igraph)

gffFile <- "ENSG00000214548.gff"
gffFile <- "MSTRG.64414.gff"
gffFile <- "/tmp/ENSG00000223972.5.gtf"
gffFile <- "ENSG00000230630.5.gtf"
#gffFile <- "ENSG00000047410.14.gtf"

gr <- import(gffFile)
gr <- gr[!(mcols(gr)$type %in% c("start_codon", "stop_codon"))]

## Create gene_id and transcript_id metadata columns:

transcript_id <- mcols(gr)$transcript_id
## Sanity check:
#stopifnot(identical(runLength(Rle(transcript_id)),
#                    runLength(Rle(mcols(gr)$name))),
#          !anyDuplicated(runValue(Rle(transcript_id))),
#          !anyDuplicated(runValue(Rle(mcols(gr)$name))))

txdb <- makeTxDbFromGRanges(gr)



SG2igraph <- function(geneID, sg, edges_by_gene) {

  nodes = sgnodes(sg[geneID])
  g1 = edges_by_gene[[geneID]]

  g1.df <- cbind.data.frame( from=g1$from, to=g1$to, g1[,4:5])
  g1.ncol = ncol(g1.df)
  tx_ids = unique(unlist(g1.df$tx_id))
  for (tx_id in tx_ids) {
    g1.df = cbind.data.frame(g1.df, unlist(lapply(lapply(g1.df$tx_id, '%in%', tx_id) , any)))
  }
  colnames(g1.df)[(g1.ncol+1):(g1.ncol+length(tx_ids))] = tx_ids

  if (g1.df$strand[1] == '+') {
    node1 <- data.frame(oldid = g1.df$from, newcoord = g1.df$start)
    node2 <- data.frame(oldid = g1.df$to, newcoord = g1.df$end+1)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(node_coord$newcoord),]
  }
  else if (g1.df$strand[1] == '-') {
    start = g1.df$start
    end = g1.df$end
    g1.df$start = end
    g1.df$end = start
    node1 <- data.frame(oldid = g1.df$from, newcoord = g1.df$start+1)
    node2 <- data.frame(oldid = g1.df$to, newcoord = g1.df$end)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(-node_coord$newcoord),]
  }
  nodes.df = data.frame( ID=c('R', node_coord$newcoord, 'L'))

  rownames(node_coord) = node_coord$oldid
  g1.df$from = node_coord[g1.df$from,"newcoord"]
  g1.df$to = node_coord[g1.df$to,"newcoord"]
 
  write.table(g1.df, paste0(geneID, ".edgetable.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  drops <- c("seqnames","strand", "tx_id")
  g1.df = g1.df[ , !(names(g1.df) %in% drops)]
  g <- graph.data.frame(g1.df, directed=TRUE, vertices=nodes.df)

  roots = nodes.df$ID[!(nodes.df$ID %in% c(g1.df$to, "R", "L"))] 
  leaves = nodes.df$ID[!(nodes.df$ID %in% c(g1.df$from, "R", "L"))]
  for (root in roots) {
    g <- g %>% add_edges(c("R", root))
  }
  for (leaf in leaves) {
    g <- g %>% add_edges(c(leaf, "L"))
  }
  return (g)
}


  sg <- SplicingGraphs(txdb,  min.ntx=1)
  edges_by_gene <- sgedgesByGene(sg)

  geneIDs <- names(edges_by_gene)
  for (geneID in geneIDs) {
    sgigraph = SG2igraph(geneID, sg, edges_by_gene)
    write_graph(sgigraph, paste0(geneID, ".graphml"), "graphml")

    pdf(paste0(geneID, ".graphml.pdf")) 
    plot(sgigraph, vertex.size=4, edge.width=1, edge.color=as.factor(E(sgigraph)$ex_or_in))
    dev.off()
  }


