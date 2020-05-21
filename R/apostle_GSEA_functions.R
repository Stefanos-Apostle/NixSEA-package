geneid_correct <- function(RNK_file){
  dup = which(duplicated(RNK_file$GeneName) == TRUE)
  dup_gn = c()
  for (i in dup) { dup_gn = c(dup_gn, RNK_file$GeneName[i]) }
  dup_gn = unique(dup_gn)

  gn_pos_list = c()
  gn = c()
  rk = c()
  for (i in dup_gn) {
    gn_pos = grep(i, RNK_file$GeneName)
    gn_pos_list = c(gn_pos_list, gn_pos)
    gn_mean = mean(RNK_file$rank[gn_pos])
    gn = c(gn, i)
    rk = c(rk, gn_mean)
  }
  avg_df = data.frame(GeneName = gn,rank = rk)
  RNK_file = RNK_file[-gn_pos_list, ] ##remove duplicated variables
  RNK_file = rbind(RNK_file, avg_df)
  RNK_file = RNK_file[order(abs(RNK_file$rank), decreasing = TRUE),]
  return(RNK_file)
}

entrez_correct <- function(RNK_file){
  dup = which(duplicated(RNK_file$entrez_id) == TRUE)
  dup_gn = c()
  for (i in dup) { dup_gn = c(dup_gn, RNK_file$entrez_id[i]) }
  dup_gn = unique(dup_gn)

  gn_pos_list = c()
  gn = c()
  rk = c()
  for (i in dup_gn) {
    gn_pos = grep(i, RNK_file$entrez_id)
    gn_pos_list = c(gn_pos_list, gn_pos)
    gn_mean = mean(RNK_file$rank[gn_pos])
    gn = c(gn, i)
    rk = c(rk, gn_mean)
  }
  avg_df = data.frame(entrez_id = gn,rank = rk)
  RNK_file = RNK_file[-gn_pos_list, ] ##remove duplicated variables
  RNK_file = rbind(RNK_file, avg_df)
  RNK_file = RNK_file[order(abs(RNK_file$rank), decreasing = TRUE),]
  return(RNK_file)
}


RNK_make <- function(deseq2_res) {
  rank_stat <- -log10(deseq2_res$padj)/sign(deseq2_res$log2FoldChange)
  RNK_file <- data.frame(entrez_id = deseq2_res$entrez, rank = rank_stat)
  length(RNK_file$entrez_id)
  pseudo <- grep("Gm", deseq2_res$mgi_symbol)
  entrez_na <- which(is.na(RNK_file$entrez_id) == TRUE)
  rank_na <- which(is.na(RNK_file$rank) == TRUE)
  del <- c(pseudo, entrez_na, rank_na)
  del <- unique(del)
  RNK_file <- RNK_file[-del,]
  RNK_file <- RNK_file[order(abs(RNK_file$rank), decreasing= TRUE),]
  RNK_file <- entrez_correct(RNK_file) ##averages duplicated entrez_genes
  fgsea_RNK <- RNK_file$rank
  names(fgsea_RNK) <- RNK_file$entrez_id
  return(fgsea_RNK)
}


msigdb_sets <- setClass("msigdb_sets", slots = c(gs = "tbl_df",
                                                 gs_names = "character",
                                                 gs_ids = "character"))

keyword_geneset <- function(database, keyword) {
  gs <- database[grep(keyword, database$gs_name),]
  gs_names <- unique(gs$gs_name)
  gs_ids <- unique(gs$gs_id)
  sets <- msigdb_sets()
  sets@gs <- gs
  sets@gs_names <- gs_names
  sets@gs_ids <- gs_ids
  return(sets)
}

geneset_list <- function(set){ #set is S4 class object created from output of function msigdb_keyword()
  i = 1
  f = c()
  for (i in 1:length(set@gs_ids)) {
    f[[set@gs_names[i]]] = set@gs$entrez_gene[which(set@gs$gs_id == set@gs_ids[i])]
    i = i+1
  }
  return(f)
}

enrichment_figures <- function(pathway, stats, sig_gs_names, title,
                               gseaParam=1,
                               ticksSize=0.2) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))

  fgsea_df = data.frame()
  i = 1
  while (i <= length(pathway)) {
    pathway_reord <- unname(as.vector(na.omit(match(pathway[[i]], names(statsAdj)))))
    pathway_reord <- sort(pathway_reord)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway_reord,
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway_reord - 1, pathway_reord))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0), z = sig_gs_names[i])
    fgsea_df <- rbind(fgsea_df, toPlot)
    diff <- (max(tops) - min(bottoms)) / 8

    i = i + 1
  }

  x=y=NULL
  g <- ggplot(fgsea_df, aes(x=x, y=y, group = z, color = z)) +
    geom_point( size=0.1) +
    geom_hline(yintercept=0, colour="black") +
    geom_line() + theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_viridis_d(name = "Gene Set") +
    theme_classic() +
    xlab(NULL) +
    ylab("Enrichment Score") +
    ggtitle(title)


  h <- ggplot(data = fgsea_df, aes(x = x, y = z, color = z, group = z)) +
    geom_point(pch = "|", cex = 3) +
    ylab("") +
    xlab("Gene Rank") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_viridis_d() +
    theme_classic() +
    theme(axis.text.y = element_blank()) +
    theme(legend.position = "none")

  egg::ggarrange(g, h, heights = c(0.75, 0.25))
}


enrichment_analysis <- function(geneset_list, fgsea_RNK, msigdb_sets, figure_header) {
  res <- fgsea(geneset_list, fgsea_RNK, nperm = 10000)
  sig_genesets <- which(res$padj < 0.05)
  if (length(sig_genesets) > 0) {
    sig_gs_names <- msigdb_sets@gs_names[sig_genesets]
    geneset_sig <- geneset_list[sig_genesets]
    figure1 <- enrichment_figures(geneset_sig, fgsea_RNK, sig_gs_names, title = figure_header)
    figure1
  }
  else {
    print("No significantly enriched gene sets :/")
  }
}



LE_heatmap <- function(geneset_list, fgsea_RNK, deseq_res, heatmap_header) {
  res <- fgsea(geneset_list, fgsea_RNK, nperm = 10000)
  sig_genesets <- which(res$padj < 0.05)
  leading_edge <- deseq_res$hgnc_symbol[which(deseq_res$entrez %in% unique(unlist(res[sig_genesets]$leadingEdge)))] #change hgnc_symbol based on organism
  hm <- as.matrix(rld_df[which(rownames(rld_df) %in%
                                 leading_edge),])
  plot.new()
  figure2 <- heatmap.2(x = hm, scale="row",
                       dendrogram= "column",
                       trace = "none",
                       col = colorRampPalette(c("darkblue", "lightblue","white", "orange", "darkred"))(n = 99),
                       labCol = colnames(rld_df),
                       #ColSideColors = col_colors,
                       ylab = "Gene",
                       xlab = "Biological Sample",
                       main = heatmap_header)
  figure2
}

































