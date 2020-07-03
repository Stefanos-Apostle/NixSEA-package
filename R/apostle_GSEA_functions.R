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


RNK_make <- function (deseq2_res, countdata, coldata,  conditions, stat="Signal2Noise") {
  if (stat == "Signal2Noise"){
    cond1 <- coldata$name[which(coldata$condition == conditions[1])]
    cond2 <- coldata$name[which(coldata$condition == conditions[2])]
    condition1_counts <- countdata[,which(colnames(placenta_counts) %in% cond1)]
    condition2_counts <- countdata[,which(colnames(placenta_counts) %in% cond2)]

    condition1_counts$sd <- rowSds(condition1_counts)
    condition2_counts$sd <- rowSds(condition2_counts)

    condition1_counts$mean <- rowMeans(condition1_counts[,1:length(cond1)])
    condition2_counts$mean <- rowMeans(condition2_counts[,1:length(cond2)])

    rank_stat <- (condition1_counts$mean - condition2_counts$mean)/(condition1_counts$sd + condition2_counts$sd) ## (mu(a)-mu(b))/(sd(a)+sd(b))

  }else if (stat == "Pl2FC") {
    rank_stat <- -log10(deseq2_res$padj)/sign(deseq2_res$log2FoldChange)
  }else{
    stop("Rank stat must be 'Signal2Noise' or 'Pl2FC'")
  }

  RNK_file <- data.frame(entrez_id = deseq2_res$entrez, rank = rank_stat)
  length(RNK_file$entrez_id)
  if ("hgnc_symbol" %in% colnames(deseq2_res)) {
    symbol <- deseq2_res$hgnc_symbol
  }
  if ("mgi_symbol" %in% colnames(deseq2_res)) {
    symbol <- deseq2_res$mgi_symbol
  }
  pseudo <- grep("Gm", symbol)
  entrez_na <- which(is.na(RNK_file$entrez_id) == TRUE)
  rank_na <- which(is.na(RNK_file$rank) == TRUE)
  del <- c(pseudo, entrez_na, rank_na)
  del <- unique(del)
  RNK_file <- RNK_file[-del, ]
  RNK_file <- RNK_file[order(abs(RNK_file$rank), decreasing = TRUE),
                       ]
  RNK_file <- entrez_correct(RNK_file)
  fgsea_RNK <- RNK_file$rank
  names(fgsea_RNK) <- RNK_file$entrez_id
  return(fgsea_RNK)
}


msigdb_sets <- setClass("msigdb_sets", slots = c(gs = "tbl_df",
                                                 gs_names = "character",
                                                 gs_ids = "character"))

##Updated function to support multiple keywords
keyword_geneset <- function (database, keyword)
{
  if (length(keyword) > 1) {
    x <- c()
    for (i in keyword) {
      y <- grep(toupper(i), database$gs_name)
      x <- c(x,y)
      x <- unique(x)
    }
  }else{
    x <- grep(keyword, database$gs_name)
  }
  gs <- database[x, ]
  gs_names <- unique(gs$gs_name)
  gs_ids <- unique(gs$gs_id)
  sets <- msigdb_sets()
  sets@gs <- gs
  sets@gs_names <- gs_names
  sets@gs_ids <- gs_ids
  return(sets)
}

#keyword_geneset <- function(database, keyword) {
#  gs <- database[grep(keyword, database$gs_name),]
#  gs_names <- unique(gs$gs_name)
#  gs_ids <- unique(gs$gs_id)
#  sets <- msigdb_sets()
#  sets@gs <- gs
#  sets@gs_names <- gs_names
#  sets@gs_ids <- gs_ids
#  return(sets)
#}

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



LE_heatmap <- function (geneset_list, fgsea_RNK, deseq_res, expression_matrix, ColumnColors = FALSE,
                        heatmap_header = NULL, col_margin = 10, row_margin = 10)
{
  res <- fgsea(geneset_list, fgsea_RNK, nperm = 10000)
  sig_genesets <- which(res$padj < 0.05)
  if ("hgnc_symbol" %in% colnames(deseq_res)) {
    symbol <- deseq_res$hgnc_symbol
  }
  if ("mgi_symbol" %in% colnames(deseq_res)) {
    symbol <- deseq_res$mgi_symbol
  }
  leading_edge <- symbol[which(deseq_res$entrez %in% unique(unlist(res[sig_genesets,
                                                                       1:length(res)]$leadingEdge)))]
  le_genes <- which(rownames(expression_matrix) %in% leading_edge)
  df <- expression_matrix[le_genes, ]
  hm <- as.matrix(df)
  plot.new()
  if (ColumnColors == TRUE) {
    i = 1
    colcolors <- as.character(dds$treatment)
    while (i <= length(levels(dds$treatment))) {
      colcolors <- replace(colcolors, grep(levels(dds$treatment)[i],
                                           colcolors), viridis(length(levels(dds$treatment)))[i])
      i = i + 1
    }
    print(colcolors)
    heatmap.2(x = hm, scale = "row", dendrogram = "column",
              trace = "none", col = colorRampPalette(c("darkblue",
                                                       "lightblue", "white", "orange", "darkred"))(n = 99),
              labCol = colnames(expression_matrix), ylab = "Gene", xlab = "Biological Sample",
              main = heatmap_header, ColSideColors = colcolors,
              margins = c(col_margin, row_margin))
    legend(x = "left", legend = as.character(levels(dds$treatment)),
           fill = viridis(length(levels(dds$treatment))))
  }else if (class(ColumnColors) == "character"){
    heatmap.2(x = hm, scale = "row", dendrogram = "column", ColSideColors = ColumnColors,
              trace = "none", col = colorRampPalette(c("darkblue",
                                                       "lightblue", "white", "orange", "darkred"))(n = 99),
              labCol = colnames(expression_matrix), ylab = "Gene", xlab = "Biological Sample",
              main = heatmap_header)

  }else {
    heatmap.2(x = hm, scale = "row", dendrogram = "column",
              trace = "none", col = colorRampPalette(c("darkblue",
                                                       "lightblue", "white", "orange", "darkred"))(n = 99),
              labCol = colnames(expression_matrix), ylab = "Gene", xlab = "Biological Sample",
              main = heatmap_header)
  }
}

## Gene set distribution plot function
geneset_dist <- function(gsea_res){
  col_fun <- function(x){if (x < 0.05){"red"}else{"black"} }
  col <- sapply(gsea_res$padj, col_fun)

  ggplot(data=gsea_res, aes(x=NES, y=padj)) +
    geom_point( color=col) +
    xlab("Normalized Enrichment Score") +
    ylab("Adjusted P value") +
    ggtitle("Gene Set Enrichment Distribution") +
    theme_classic()
}

## Visualization of top gene sets enriched in each condition
top_genesets <- function (gsea_res, top_num = 5, conditions = c()) {
  gsea_up <- gsea_res[which(gsea_res$NES > 0), ]
  gsea_up <- gsea_up[order(gsea_up$NES, decreasing = TRUE),
                     ]
  gsea_down <- gsea_res[which(gsea_res$NES < 0), ]
  gsea_down <- gsea_down[order(gsea_down$NES, decreasing = FALSE),]

  if (length(conditions) == 0) {
    up_cols <- c("Enriched in Up Condition", "P-adj","NES")
    down_cols <- c("Enriched in Down Condition", "P-adj","NES")
  }else if (length(conditions) == 2) {
    up_cols <- c(paste("Enriched in", conditions[1], sep = " "), "P-adj","NES")
    down_cols <- c(paste("Enriched in", conditions[2], sep=" "), "P-adj","NES")
  }else{
    stop("conditions must be empty or a list of length two with your GSEA condition names.")
  }

  up <- ggtexttable(x = gsea_up[1:top_num, c(1, 3, 5)], rows = NULL,
                    theme = ttheme("lRedWhite"), cols = up_cols)
  down <- ggtexttable(x = gsea_down[1:top_num, c(1, 3, 5)],
                      rows = NULL, theme = ttheme("lBlueWhite"), cols = down_cols)
  ggarrange(up, down, nrow = 2)
}

## Annotation of keywords in significantly enriched gene sets
keyword_ann <- function (gsea_res, top_keywords = 10) {
  sig_paths <- length(which(gsea_res$padj < 0.05))
  all_words <- unlist(strsplit(gsea_res$pathway[order(gsea_res$padj, decreasing = FALSE)][1:sig_paths],
                               "_"))
  un_words <- unique(all_words)
  common_words <- c("THE", "HALLMARK", "AND", "OR", "UP", "DN",
                    "GENES", "REACTOME", "KEGG", "GO", "PATHWAY", "VIA",
                    "WITH", "CELL", "BIOCARTA")
  n_rep <- c()
  short <- c()
  for (i in un_words) {
    if (nchar(i) > 2) {
      l <- length(which(all_words %in% i))
      n_rep <- c(n_rep, l)
    }
    else {
      s <- match(i, un_words)
      short <- c(short, s)
    }
  }
  ann_word <- data.frame(Keyword = un_words[-short], Number = n_rep)
  ann_word <- ann_word[order(ann_word$Number, decreasing = TRUE),
                       ]
  cm <- which(ann_word$Keyword %in% common_words)
  ann_word <- ann_word[-cm, ]
  num_gs <- top_keywords
  avg_p <- c()
  le_k <- c()
  for (i in ann_word$Keyword[1:num_gs]) {
    paths <- grep(i, gsea_res$pathway[gsea_res$padj < 0.05])
    m <- mean(gsea_res$padj[gsea_res$padj < 0.05][paths])
    avg_p <- c(avg_p, m)
    k <- length(unique(gsea_res$leadingEdge[paths]))
    le_k <- c(le_k, k)
  }
  ggplot(data = ann_word[1:num_gs, ]) + geom_point(aes(x = (Number/sig_paths) *
                                                         100, y = fct_rev(factor(Keyword, level = Keyword)), color = avg_p,
                                                       size = le_k)) + scale_color_viridis() + ylab("Top Keywords") +
    xlab("Percent of Significant Gene Sets Including Keyword") +
    labs(size = "Unique LE Genes", color = "Average padj") +
    ggtitle("Significant Gene Set Keyword Annotation") +
    theme_classic()
}

##function to average mgi symbol duplicates
mgi_dup_fix <- function(rld_df, num_samples){
  un_syms <- unique(rld_df$mgi_symbol)[-which(is.na(rld_df$mgi_symbol) == TRUE)]
  revised_df <- data.frame()
  for (i in un_syms){
    pos <- which(rld_df$mgi_symbol %in% i)
    x <- data.frame()
    for (j in pos){
      if (nrow(x) == 0) {
        x <- rbind(x,rld_df[j, 1:num_samples])
      }else{
        x <- x + rld_df[j, 1:num_samples]
      }
    }
    avg_x <- x/length(pos)
    if (nrow(revised_df) == 0){
      revised_df <- as.data.frame(avg_x)
    }else{
      revised_df <- rbind(revised_df, avg_x)
    }
  }
  rownames(revised_df) <- un_syms
  return(revised_df)
}





























