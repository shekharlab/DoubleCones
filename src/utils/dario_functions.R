# Dario's functions!

custom_median <- function(x) {
  x <- sort(x, na.last = NA)  # Sort the values and remove NAs
  n <- length(x)
  
  if (n %% 2 == 1) {
    # If odd, return the middle element
    return(x[(n + 1) / 2])
  } else {
    # If even, return the lower of the two middle elements
    return(x[n / 2])
  }
}

LegendLowerRight = function(){
  theme(
    legend.position = c(0.95, 0.05),  # Place the legend inside at the lower right
    legend.justification = c(1, 0)  # Adjust the anchor point of the legend
  )  
}

#' Install a Package from CRAN or Bioconductor
#'
#' This function attempts to install a given package from CRAN. If the installation
#' from CRAN fails, it will attempt to install the package from Bioconductor using `BiocManager`.
#' If the package is already installed, the function will skip the installation.
#'
#' @param pkg A character string specifying the name of the package to install.
#'
#' @return The function does not return a value. It attempts to install the specified package.
#' It will print messages to indicate the success or failure of the installation attempts.
#'
#' @details
#' The function first checks if the package is already installed using `requireNamespace()`. If the package
#' is not installed, it first tries to install it from CRAN using `install.packages()`. If the installation
#' from CRAN fails (e.g., the package is only available on Bioconductor), the function will check if
#' `BiocManager` is installed. If not, it installs `BiocManager` and then attempts to install the package
#' from Bioconductor.
#'
#' @examples
#' \dontrun{
#' # Install ggplot2 from CRAN or Bioconductor if it's not already installed:
#' install_package("ggplot2")
#'
#' # Install a Bioconductor package (e.g., "GenomicFeatures"):
#' install_package("GenomicFeatures")
#' }
#'
#' @seealso \code{\link{install.packages}}, \code{\link{BiocManager::install}}
#'
#' @export
install_package <- function(pkg) {
  # Check if the package is already installed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Try to install from CRAN
    tryCatch({
      message(paste("Attempting to install", pkg, "from CRAN..."))
      install.packages(pkg)
    }, error = function(e) {
      message(paste("Failed to install", pkg, "from CRAN. Trying BiocManager..."))
      
      # Check if BiocManager is installed
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      
      # Try to install from Bioconductor
      tryCatch({
        BiocManager::install(pkg)
      }, error = function(e) {
        message(paste("Failed to install", pkg, "from Bioconductor."))
      })
    })
  } else {
    message(paste("Package", pkg, "is already installed."))
  }
}


DrawPhylogeneticTreeFromNewick = function(newick, key, as_ylab = TRUE, layout = 'rectangular', species.prune = NULL, rootedge = 50, interval = 50){
  
  # Convert to common name
  tree = ape::read.tree(newick)
  tree$tip.label = key$CommonName[match(gsub("_", " ", tree$tip.label), key$LatinName)]
  # dendro = phylogram::as.dendrogram.phylo(tree)
  # sorted.dendro = reorder(dendro, 1:20)
  
  # Prune species if necessary
  # pruned.dendro = dendextend::prune(sorted.dendro, species.prune)
  
  # phylo_obj = ape::as.phylo(as.hclust(pruned.dendro))
  phylo_obj = tree
  breaks = seq(0, max(phylo_obj$edge.length)+interval, by = interval)
  max_time = max(phylo_obj$edge.length)
  shift = max_time-breaks[length(breaks)-1]
  
  ggtree(phylo_obj, layout = layout) + 
    theme_tree2() + 
    geom_tiplab(as_ylab=as_ylab) + 
    geom_rootedge(rootedge = rootedge) +
    scale_x_continuous(name = "Time (millions of years ago)",
                       breaks = seq(shift-interval, max_time, by = 50),  # Breaks at 50 MYA intervals
                       labels = rev(breaks))+  # Adjust labels so 0 MYA is present
    theme(axis.text = element_text(color = 'black'))
}

rotate_coords <- function(coords, angle, radians = TRUE) {
  # Convert the angle from degrees to radians
  if(!radians) angle_radians <- angle * pi / 180 else angle_radians = angle
  
  # Create the rotation matrix
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                              sin(angle_radians), cos(angle_radians)), 
                            nrow = 2, byrow = TRUE)
  
  # Original coordinates
  # coords <- cbind(x, y)
  
  # Rotate the coordinates
  rotated_coords <- t(rotation_matrix %*% t(coords))
  
  # Return the rotated coordinates
  return(data.frame(rotated_coords[, 1], rotated_coords[, 2]) %>% setNames(colnames(coords)))
}

UMIplot <- function(counts, cutoff = 200, title = NULL){
  count_data=data.frame(rank=1:length(counts), 
                        count=rev(sort(counts)))
  cutoff.rank=nrow(count_data[count_data$count>=cutoff,])
  # plot=plot(log10(count+1) ~ rank, data=count_data, log="x", type = "S", main=title)
  # abline(v=cutoff.rank, lty=2)
  ggplot(count_data, aes(x = rank, y = log10(count + 1))) +
    geom_step() +  # equivalent to type = "S"
    scale_x_log10() +  # equivalent to log = "x"
    geom_vline(xintercept = cutoff.rank, linetype = 'dashed') +
    labs(title = title, x = "Rank", y = "log10(Count + 1)") +
    theme_bw()
}

MakeMutuallyExclusive = function(list){
  all = unlist(list)
  duplicates = getDuplicates(all)
  lapply(list, function(x) x[!x %in% duplicates])
}

setRowNames <- function(object, names) {
  rownames(object) <- names
  return(object)
}

# Function to compute softmax
softmax <- function(x) {
  exp_x <- exp(x)        # Exponentiate each element
  return(exp_x / sum(exp_x))  # Normalize by the sum of exponentiated elements
}

GetTrainingSet = function(object, group.by, proportion = 0.6, max.size = 100, seed = 12345){
  
  Idents(object) = group.by
  clusters = factor(object@meta.data[[group.by]])
  training_cells = unlist(sapply(seq_along(levels(clusters)), function(index) 
    WhichCells(object, idents = levels(clusters)[[index]], 
               downsample = floor(table(clusters)[[index]])*proportion)))
  object_train = object[,training_cells]
  object_train = DownsampleSeurat(object_train, group.by = group.by, size = max.size, seed = seed)
  return(object_train)
}

PieChart = function(df, label = FALSE, label.func = geom_text_repel){
  ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
    geom_col(color = "black") +
    {if(label) label.func(aes(label = Var1), position = position_stack(vjust = 0.5))} +
    coord_polar(theta = "y") +
    # scale_fill_brewer() +
    labs(fill = "Cluster") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.background = element_rect(fill = "white")) 
}

FindGenes = function(object, prefix){
  rownames(object)[startsWith(rownames(object), prefix)]
}

DimPlotLabeled = function(object, group.by = 'seurat_clusters', ...){
  DimPlot(object, group.by = group.by, label = TRUE, ...) + NoLegend()
}

convert_values <- function(data_vector, key_table) {
  # Ensure the key_table has the correct columns
  if (!all(c("old.names", "new.names") %in% names(key_table))) {
    stop("key_table must have 'old.names' and 'new.names' columns")
  }
  
  # Perform the conversion
  converted_vector <- key_table$new.names[match(data_vector, key_table$old.names)]
  
  # Return the converted vector
  return(converted_vector)
}

MoreThanXCells = function(object, annotate.by = 'seurat_clusters', group.by = NULL, more.than = 0){
  tabulation = table(object@meta.data[[annotate.by]], object@meta.data[[group.by]])
  print(tabulation)
  names(which(apply(as.data.frame.matrix(tabulation), 1, function(x) all(x > more.than))))
}

exchange_factor_level <- function(factor_var, old_level, new_level) {
  # Check if the input is a factor
  if (!is.factor(factor_var)) {
    stop("The input variable must be a factor.")
  }
  
  # Check if the old level exists in the factor
  if (!(old_level %in% levels(factor_var))) {
    stop("The old level does not exist in the factor.")
  }
  
  # Check if the new level already exists in the factor
  if (new_level %in% levels(factor_var)) {
    stop("The new level already exists in the factor.")
  }
  
  # Replace the old level with the new level
  levels(factor_var)[levels(factor_var) == old_level] <- new_level
  
  return(factor_var)
}

geneDotPlotStack = function(object, genes){
  StackedPlots(lapply(genes, geneDotPlotFast))
}

PlainYAxis = function(){
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())
}

PlainXAxis = function(){
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())
}

get_legend <- function(my_plot) {
  plot_grob <- ggplotGrob(my_plot)
  legend <- plot_grob$grobs[[which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

add_nulls_between <- function(lst) {
  # Use lapply to insert NULLs
  interleaved <- lapply(lst, function(x) list(x, NULL))
  # Flatten the list using do.call and c
  flattened <- do.call(c, interleaved)
  # Remove the last NULL added by the function
  return(flattened[-length(flattened)])
}

smartReadRDS = function(filepath){
  objectname = basename(gsub('.rds', '', filepath))
  if(!exists(objectname)){
    readRDS(filepath)
  } else {
    message(objectname, ' already loaded')
  }
}

VlnPlot2 = function(object, features = NULL, split.by = NULL, group.by = NULL, fill.by = 'feature', cols = NULL, idents = NULL, stack = FALSE, pt.size = 0, combine = TRUE){
  meta = Metadata(object, split.by, group.by)
  VlnPlot(object, features = features, split.by = split.by, group.by = group.by, stack = stack, 
          idents = idents, fill.by = fill.by, cols = cols[match(meta[[group.by]], names(cols))], pt.size = pt.size, combine = combine)
}

RunDESeq2 = function(object, ident.1, ident.2){
  
  colData = object@meta.data
  colData$condition = factor(colData$condition, levels = c(ident.2, ident.1))
  dds <- DESeqDataSetFromMatrix(countData = object@assays$RNA@counts,
                                colData = colData,
                                design= ~ condition)
  dds <- DESeq(dds)
  result_name = resultsNames(dds)[2]
  res <- results(dds, name=result_name)
  result = data.frame(gene = rownames(res), 
                      cluster = str_split_fixed(result_name, '_', 4)[,2], 
                      avgExpr = res$baseMean, 
                      avg_log2FC = res$log2FoldChange, 
                      stat = res$stat,
                      p_val = res$pvalue,
                      p_val_adj = res$padj) %>% arrange(p_val_adj, avg_log2FC)
  
  # sorted.res = res %>% 
  # as.data.frame %>% 
  # setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out')) %>% 
  # arrange(padj, log2FoldChange)
  
  result
}

DendroOrder2 = function(object, group.by = 'seurat_clusters', ...){
  dendro = CorrelationHeatmap(object, group.by = group.by, return.dendrogram = TRUE, ...)
  object$dendro.order = factor(object@meta.data[[group.by]], levels = labels(dendro))
  object@tools$BuildClusterTree = ape::as.phylo(dendro)
  # PlotClusterTree(object)
  return(object)
}

ConvertGeneSymbolsFromMatrix = function(matrix, list){
  subset_list = list[list %in% colnames(matrix)]
  subset_matrix = matrix[,match(subset_list, colnames(matrix))]
  # column_vector = as.numeric(colnames(subset_matrix) %in% list)
  new_genes = unlist(apply(subset_matrix, 2, function(x) names(which(x == 1))))
  return(new_genes)
}

CelltypeProportionBarplotGenotype = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'genotype'){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(seurat, 'sample_core', 'tissue', 'treatment', 'genotype')
  prop_data = cbind(melted, metadata[match(melted$Sample, metadata$sample_core),-1])
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ genotype) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

CelltypeProportionBarplotTreatment = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'treatment'){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(seurat, 'sample_core', 'tissue', 'treatment', 'genotype')
  prop_data = cbind(melted, metadata[match(melted$Sample, metadata$sample_core),-1])
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ treatment) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .2))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

CelltypeProportionBarplotTissue = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'tissue'){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(seurat, 'sample_core', 'tissue', 'treatment', 'genotype')
  prop_data = cbind(melted, metadata[match(melted$Sample, metadata$sample_core),-1])
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ tissue) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

getDuplicates = function(vector){
  unique(vector[duplicated(vector)])
}

SeuratHeatmap = function(object, features, group.by = 'seurat_clusters', color = 'red', max.z.score = 2, assay = 'RNA', rotate = FALSE, annotate.by = NULL, anno_cols = NULL, ...){
  
  # scaled.expr = t(scale(t(expm1(LogAvgExpr(object, assay = assay, features = features)))))
  
  # Using scale.data slot because this is at the cell level, not averaged
  scaled.expr = as.data.frame(object@assays[[assay]]@scale.data)
  colData = object@meta.data
  genes = rownames(scaled.expr)
  colData = colData[order(colData[[group.by]]),]
  scaled.expr = scaled.expr[match(features, genes),rownames(colData)]
  
  if(!is.null(annotate.by)){
    metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by)
    annotation = metadata[match(object$species, metadata[[group.by]]), annotate.by]
    if(is.null(anno_cols)){
      ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left")
    } else{
      ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left", col = list(anno = annotation_cols))
    }
  } else {
    ha = NULL
  }
  
  if(rotate) scaled.expr = t(scaled.expr)
  
  Heatmap(as.matrix(scaled.expr),
          name = "mat",
          col = colorRamp2(c(-max.z.score, 0, max.z.score), c("white", "white", color)),
          cluster_rows = FALSE,
          cluster_columns = FALSE, 
          # row_split = factor(row_split, levels = unique(row_split)), 
          border = TRUE, 
          heatmap_legend_param = list(title = ""),
          # rect_gp = gpar(col = "black", lwd = 1),
          # left_annotation = rowAnnotation(#row = anno_textbox(gene.key$row_annotation, gene.key$row_annotation), 
          #                     row = anno_block(labels = unique(gene.key$row_annotation), 
          #                                labels_gp = gpar(col = "white", fontsize = 10))),
          top_annotation = ha, 
          ...)
  
  # Add horizontal bars separating rows
  # breaks = (sapply(unique(types), function(x) which(types == x)[1])-1)[-1]
  # for(index in breaks) {
  #   for(slice in seq_along(unique(row_split))){
  #     x_coord = index/length(types)
  #     
  #     # Vertical
  #     decorate_heatmap_body("mat", row_slice = slice, {
  #       grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
  #     })
  #   }
  # }
  
}

volcanoPlot <- function(table, fdr_cutoff = 0.05, fc_cutoff = 0.25, max_fdr = 1e-100, max_fc = 2, labels = FALSE, ...){
  # max_fc = max(abs(table$avg_log2FC))
  table$p_val_adj[table$p_val_adj < max_fdr] = max_fdr
  table$avg_log2FC[table$avg_log2FC < -max_fc] = -max_fc
  table$avg_log2FC[table$avg_log2FC > max_fc] = max_fc
  data_signif = subset(table, p_val_adj < fdr_cutoff & abs(avg_log2FC) > fc_cutoff)
  data_ns = subset(table, p_val_adj >= fdr_cutoff | abs(avg_log2FC) <= fc_cutoff)
  plot=ggplot(table, aes(y=(-log10(p_val_adj)), x=avg_log2FC))+
    geom_point()+
    ylab("-log10 FDR")+
    xlab("log2 fold change")+
    geom_point(data = data_signif, color = "red")+
    geom_point(data = data_ns, color = "black")+
    {if(labels) geom_text_repel(data = data_signif, label = data_signif$gene, ...)}+
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 'dashed')+
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 'dashed')+
    annotation_custom(
      grob = textGrob(nrow(subset(data_signif, avg_log2FC < 0)), x = 0.1, y = 0.2, hjust = 1, vjust = 1, gp = gpar(col = "red")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    annotation_custom(
      grob = textGrob(nrow(subset(data_signif, avg_log2FC > 0)), x = 0.9, y = 0.2, hjust = 1, vjust = 1, gp = gpar(col = "red")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    # annotate('text', x = 0.2, y = 0.8, label = nrow(subset(data_signif, avg_log2FC < 0)))+
    scale_x_continuous(limits = c(-max_fc, max_fc))+
    theme_cowplot()
  return(plot)
}

FindMarkersFast = function(object, group.by = 'seurat_clusters', ident.1 = NULL, ident.2 = NULL, p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 0.25){
  output = wilcoxauc(object, group_by = group.by, groups_use = c(ident.1, ident.2)) %>% setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
  if(is.factor(object@meta.data[[group.by]])) output$cluster = factor(output$cluster, levels = levels(object@meta.data[[group.by]]))
  output %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff & cluster == ident.1) %>% arrange(cluster, p_val_adj, desc(abs(avg_log2FC)))
}

FindAllMarkersFast = function(object, group.by = 'seurat_clusters', p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 0.25){
  output = wilcoxauc(object, group_by = group.by) %>% setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
  if(is.factor(object@meta.data[[group.by]])) output$cluster = factor(output$cluster, levels = levels(object@meta.data[[group.by]]))
  output %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff) %>% arrange(cluster, p_val_adj, -avg_log2FC)
}

JSHeatmap2 = function(list1, list2, xlab = NULL, ylab = NULL, ...){
  args = list(...)
  JSHeatmap(JSMatrix(table((list1), as.character(list2))), ...)+
    {if(is.null(args$title)) ggtitle(paste0('ARI = ', signif(adj.rand.index(list1, list2), 2)))}+ 
    xlab(xlab)+
    ylab(ylab)+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'))
}

ConfusionMatrix = function(true, predicted, row.norm = TRUE, plot = TRUE, max.value = 1, col.high = "#584B9FFF", col.low = 'white', legend_name = 'Percentage\nmapping'){
  raw_table = table(true, predicted)
  if(row.norm){
    matrix = raw_table/rowSums(raw_table)
  } else {
    matrix = t(t(raw_table)/colSums(raw_table))
  }
  
  melted = reshape2::melt(matrix) %>% setNames(c('True', 'Predicted', 'Percentage'))
  plt = ggplot(melted, aes(x = Predicted, y = True))+
    geom_point(aes(colour = Percentage,  size=Percentage))+
    scale_color_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = 'grey') +
    scale_radius(legend_name, limits=c(0, max.value)) +
    # scale_size(range = c(1, max.size), limits = c(0, max.perc))+
    theme_bw() +
    RotatedAxis() +
    ylab('True')+
    xlab('Predicted')+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5))
  
  if(plot) return(plt) else return(matrix)
}

AddSpacing = function(plt.list){
  
}

ConvertGeneSymbols2 = function(object, orthology_graph){
  
  # Convert matrix into a key
  key = as.data.frame(ConvertGeneSymbolsFromMatrix(orthology_graph, rownames(object))) %>% rownames_to_column() %>% setNames(c('old.names', 'new.names'))
  
  # Save original feature names
  original.features = rownames(object)
  object@misc$orig.features = original.features
  
  # Rename features
  converted.features = convert_values(rownames(object), key)
  converted.features[is.na(converted.features)] = original.features[is.na(converted.features)]
  object = RenameFeatures(object, new.names = converted.features)
  
  return(object)
}

PlotBinaryMask = function(full_mask, mask){
  grid.grabExpr(draw(Heatmap(rbind(full_mask, mask), 
                             cluster_rows = F, 
                             cluster_columns = F, 
                             col = c("grey", "blue"), 
                             column_names_side = "bottom", 
                             row_names_side = "left", 
                             rect_gp = gpar(col = "black", lwd = 1),
                             name = "ProbExp                  "
  )))
}

GrabSamples = function(filename, column){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  return(do.call(cbind, lapply(matrix.list, function(matrix) matrix[,column])))
}

StdTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  sd = apply(simplify2array(matrix.list), 1:2, sd)
  return(sd)
}

MedianTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  median = apply(simplify2array(matrix.list), 1:2, median)
  return(median)
}

MeanTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  mean = apply(simplify2array(matrix.list), 1:2, mean)
  return(mean)
}

DEGTree = function(object, group.by = "annotated", title = NULL, mc.cores = 1, ...){
  types.use = unique(object@meta.data[[group.by]])
  ndegs = Iterate(types.use, function(i,j) nDEGs(wilcoxauc(object, group_by = group.by, groups_use = c(i,j)), presto = TRUE, ...), mc.cores = mc.cores)
  
  dist = matrix(, nrow = length(types.use), ncol = length(types.use), dimnames = list(types.use, types.use))
  dist[lower.tri(dist)] = unlist(ndegs)
  dist[upper.tri(dist)] = t(dist)[upper.tri(t(dist))]
  dist[is.na(dist)] = 0
  
  Heatmap(dist, 
          col = colorRamp2(c(0, max(dist, na.rm = TRUE)), c("white", "red")),
          cluster_rows = hclust(as.dist(dist)), 
          cluster_columns = hclust(as.dist(dist)), 
          clustering_method_columns = "manhattan",
          clustering_method_rows = "manhattan",
          show_row_names = TRUE, 
          show_column_names = TRUE, 
          row_dend_side = "left", 
          show_column_dend = TRUE, 
          column_title = title,
          border_gp = gpar(col = "black", lty = 1),
          # column_title = title,
          cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.0f", dist[i, j]), x, y, gp = gpar(fontsize = 10))},
          heatmap_legend_param = list(title = "# DEGs"),
          top_annotation = NULL)
}

Iterate = function(list, FUN, return.matrix = FALSE, ...){
  outlist = lapply(seq_along(list[1:(length(list)-1)]), function(i) {
    lapply((i+1):length(list), function(j){
      FUN(list[[i]], list[[j]], ...)
    })
  })
  
  if(return.matrix) {
    matrix = matrix(, nrow = length(list), ncol = length(list))
    matrix[lower.tri(matrix)] = unlist(outlist)
    rownames(matrix) = names(list)
    colnames(matrix) = names(list)
    return(matrix) 
  } else {
    return(outlist)
  }
}

RemoveYAxis = function(plt){
  return(plt + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))
}

SeuratBootstrap = function(object, group.by = "seurat_clusters"){
  clusters = unique(object@meta.data[[group.by]])
  
  bootstrap_cols = sapply(clusters, function(cluster) {
    this_cols = which(object@meta.data[[group.by]] == cluster)
    return(sample(this_cols, length(this_cols), replace = TRUE))
  }) %>% unlist
  
  # Doesn't subset unless it's shorter than original object
  boot_object = suppressWarnings(object[,bootstrap_cols[1:length(bootstrap_cols)-1]])
  
  boot_object = RenameCells(boot_object, new.names = make.unique(Cells(boot_object)))
  
  return(boot_object)
}

leaveOneOutSupport = function(matrix, ref.tree){
  
  loo.trees.data = mclapply(seq_len(ncol(matrix)), function(col) {
    phydat = phyDat(matrix[,-col], type="USER", levels=c("0", "1"), compress = FALSE)
    tree = pratchet(phydat, trace=0, minit=100)
    list(tree = tree, parsimony = parsimony(tree, phydat))
  }, mc.cores = 16)
  
  message("Generated ", length(loo.trees.data), " LOO trees!")
  
  loo.trees = as.multiPhylo(lapply(loo.trees.data, function(x) x$tree)) %>% setNames(rownames(matrix))
  treeBS = plotBS(ref.tree, loo.trees, main="LOO support", type = "unrooted")
  ref.tree$loo.support = treeBS$node.label
  return(ref.tree)
}

FindGini = function(node){
  # TF switches at node
  tfs = gained.summary[[as.character(node)]]
  
  # Descendants
  descendants = treeRooted$tip.label[Descendants(treeRooted, as.numeric(node))[[1]] ]
  
  lapply(tfs, function(tf){
    
    data = data.table(label = orthotypes %in% descendants, expr = as.character(tf.expr.mask[tf,]))
    data.frame(tf = tf, 
               node = node, 
               gini = as.numeric(gini_impurities(data)[2,3]))
  })
}

nMonophyletic = function(...){
  length(which(FindMonophyletic(...)))
}

# For each TF, test whether it's marked OTs are monophyletic!
FindMonophyletic = function(datExpr, tree, exclude.single = TRUE){
  
  # Whether to exclude TF events specific to a single orthotype
  min = ifelse(exclude.single, 1, 0)
  max = ifelse(exclude.single, ncol(datExpr)-1, ncol(datExpr))
  
  monophyletic = sapply(rownames(datExpr), function(gene){
    orthotypes.marked = names(which(datExpr[gene,] > 0))
    ifelse(length(orthotypes.marked) > min & length(orthotypes.marked) < max, is.monophyletic(phy = tree, tips = orthotypes.marked), NA)
  })
  return(monophyletic)
}

PrettyTree = function(tree, gene, cols = c("0" = "black", "1" = "red", "0.5" = "gold")){
  # TitlePlot(
  ggtree(as.treedata(tree), layout = "equal_angle", aes(color = eval(parse(text = gene)))) + 
    geom_tippoint(aes(color=eval(parse(text = gene))), show.legend = TRUE)+
    geom_nodepoint(aes(color=eval(parse(text = gene))), show.legend = TRUE)+
    # scale_color_gradient2(low = "black", high ="red", mid = "yellow", midpoint = 0.5)+
    scale_color_manual(name = gene, values = cols)+
    # geom_tiplab(geom = "text")+
    geom_text_repel(aes(label=name), max.overlaps = Inf, show.legend = FALSE, seed = 42)+
    # geom_treescale()+
    theme(legend.position = "none")
  # gene)
}

ShuffleRows = function(matrix){
  shuffled = t(apply(matrix, 1, function(row) sample(row)))
  colnames(shuffled) = colnames(matrix)
  return(shuffled)
}

MaxParsimonyTree = function(datExpr, minit){
  phydat = phyDat(t(datExpr), type="USER", levels=c("0", "1"), compress = FALSE)
  treeBinary  <- pratchet(phydat, trace = 0, minit=minit)
  treeBinary  <- acctran(treeBinary, phydat)
  treeBinary  <- di2multi(treeBinary, tol = 1e-8)
  if(inherits(treeBinary, "multiPhylo")){
    treeBinary <- unique(treeBinary)
  }
  
  # treeBinary$node.label = (length(treeBinary$tip.label)+1):((length(treeBinary$tip.label)+treeBinary$Nnode))
  # treeBinary <- root(phy = treeBinary, node = 71)
  # plot(treeBinary)
  
  # message("Parsimony score: ", parsimony(treeBinary, phydat))
  return(list(tree = treeBinary, 
              parsimony = parsimony(treeBinary, phydat), 
              monophyletic = FindMonophyletic(datExpr, treeBinary)))
}

BinarizeExpression = function(matrix, return.vector = TRUE, plot = FALSE, seed = 42, replace.na = 0, noise = FALSE, verbose = FALSE, alpha = 0.05){
  
  library(MASS)
  library(mixtools)
  set.seed(seed)
  
  # Remove rows with zero variance, replace other NAs with replace.na
  matrix = matrix %>% replace(is.na(.), replace.na)
  input = matrix[apply(matrix, 1, var) != 0, ]
  
  vector_encoding = tryCatch({
    
    # Fit a single Gaussian distribution
    single_gaussian <- fitdistr(unlist(input), "normal")
    log_likelihood_single <- as.numeric(logLik(single_gaussian))
    
    # Fit a mixture of two Gaussian distributions
    model = if(noise) mixtools::normalmixEM(unlist(input) + rnorm(length(unlist(input)), mean = 0, sd = 0.1)) else mixtools::normalmixEM(unlist(input), k = 2)
    if(plot) plot(model, density=TRUE, breaks = 30)
    log_likelihood_mixture <- model$loglik
    
    # Perform the Likelihood Ratio Test
    LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
    p_value <- pchisq(LRT_statistic, df = 3, lower.tail = FALSE)
    if(verbose) message('LRT p-value: ', p_value)
    
    if(p_value < alpha){
      mask = matrix(model$posterior[,which.max(model$mu)], 
                    # ifelse(model$posterior[,which.max(model$mu)] > 0.8, 1, # 1 if confidently expressed, 0.5 if not confident, and 0 if confidently unexpressed
                    # ifelse(model$posterior[,which.max(model$mu)] > 0.2, 0.5, 0)),
                    nrow = nrow(input), 
                    ncol = ncol(input), 
                    dimnames = dimnames(input))
      vector_encoding = apply(mask, 2, mean)
      if(return.vector) vector_encoding else mask
    } else {
      NULL
    }
    
  }, error = function(e) {
    print(e)
  })
  vector_encoding
}

UpperCase_genes = function(object, integration = FALSE){
  
  object@misc$orig.features = rownames(object@assays$RNA@counts)
  
  rownames(object@assays$RNA@counts) <- make.unique(toupper(rownames(object@assays$RNA@counts)))
  rownames(object@assays$RNA@data)<- make.unique(toupper(rownames(object@assays$RNA@data)))
  rownames(object@assays$RNA@scale.data)<- make.unique(toupper(rownames(object@assays$RNA@scale.data)))
  object@assays$RNA@var.features <- make.unique(toupper(object@assays$RNA@var.features))
  #rownames(object@assays$RNA@meta.features)<- toupper(rownames(object@assays$RNA@meta.features))
  
  if(integration){
    rownames(object@assays$integrated@counts)<- make.unique(toupper(rownames(object@assays$integrated@counts)))
    rownames(object@assays$integrated@data)<- make.unique(toupper(rownames(object@assays$integrated@data)))
    rownames(object@assays$integrated@scale.data)<- make.unique(toupper(rownames(object@assays$integrated@scale.data)))
    object@assays$integrated@var.features <- make.unique(toupper(object@assays$integrated@var.features))
    #rownames(object@assays$integration@meta.features)<- toupper(rownames(object@assays$integration@meta.features))
  }
  
  return(object)
}

PcaContribution = function(pca, dim){
  summary = as.data.frame(pca$rotation) %>% arrange(desc(!!sym(dim)))
  contributors = summary[,dim]
  names(contributors) = rownames(summary)
  return(contributors)
}

VennDiagram = function(list1, list2, population = 17000, name1 = NULL, name2 = NULL){
  library(eulerr)
  
  summary = list(unique(list1[!list1 %in% list2]),
                 unique(list2[!list2 %in% list1]), 
                 intersect(list1,list2))
  
  # Input in the form of a named numeric vector
  input = c(
    "A" = length(unique(list1)) - length(intersect(list1,list2)),
    "B" = length(unique(list2)) - length(intersect(list1,list2)),
    "A&B" = length(intersect(list1,list2))
  )
  
  if(is.null(name1) | is.null(name2)){
    names(input) = c(deparse(substitute(list1)), deparse(substitute(list2)), paste0(deparse(substitute(list1)), "&", deparse(substitute(list2))))
  } else{
    names(input) = c(name1, name2, paste0(name1, "&", name2))
  }
  names(summary) = c(deparse(substitute(list1)), deparse(substitute(list2)), paste0(deparse(substitute(list1)), "&", deparse(substitute(list2))))
  fit <- euler(input)
  
  print(
    plot(fit, 
         quantities = TRUE,
         # fill = "transparent",
         fill = c('deepskyblue', 'orangered'),
         lty = 1,
         labels = list(font = 4))
  )
  
  # Hypergeoemtric p-value
  hyper_pval = phyper(length(intersect(list1,list2))-1, # number of white balls drawn
                      length(unique(list1)), # number of white balls
                      population-length(unique(list1)), # number of black balls
                      length(unique(list2)), # number of balls drawn
                      lower.tail=FALSE, 
                      log.p=FALSE)
  
  summary$pval = hyper_pval
  
  return(summary)
}

SAMapAlignmentHeatmap2 = function(MappingTable, type_palette, species_palette, 
                                  cols = c("white", "#fbd9d3", "#ffb09c", "red", "darkred"), 
                                  add.breaks = TRUE, show_annotation_legend = FALSE, 
                                  table_order = NULL, species.use = NULL, ...){
  
  # Cell type ordering
  species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
  types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
  idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  
  if(!is.null(species.use)){
    MappingTable = MappingTable[species %in% species.use, species %in% species.use]
    
    # re-parse
    species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
    types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
    idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  }
  
  sort_key = data.frame(species = species, types = types, idents = idents) %>% arrange(types, species)
  table_order = paste0(sort_key$species, '_', sort_key$types, '.', sort_key$idents)
  table_order = sub("\\.$", "", table_order)
  
  # table_order = as.vector(sapply(names(pr_palette), function(type) sapply(names(species_palette2), function(species) paste0(species, "_", type))))
  if(is.null(table_order)) cnames = factor(colnames(MappingTable), levels = sort(unique(colnames(MappingTable)))) else cnames = factor(colnames(MappingTable), levels = table_order)
  MappingTable = MappingTable[order(cnames), order(cnames)]
  
  # Set same species alignment scores to NA since these are never computed
  MappingTable[outer(species[order(cnames)], species[order(cnames)], FUN = function(x,y) x == y)] = NA
  
  # Type ordering
  row_km = factor(factor(types[order(cnames)], levels = names(type_palette)))
  
  print(Heatmap(MappingTable, 
                name = "Alignment\nscore", 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = HeatmapAnnotation(species = species[order(cnames)], 
                                                   `cell type` = row_km,
                                                   col = list(`cell type` = type_palette, species = species_palette), 
                                                   border = TRUE, 
                                                   show_annotation_name = TRUE, 
                                                   show_legend = FALSE), 
                left_annotation = rowAnnotation(species = species[order(cnames)], 
                                                `cell type` = row_km,
                                                col = list(`cell type` = type_palette, species = species_palette), 
                                                border = TRUE, 
                                                show_legend = show_annotation_legend,
                                                show_annotation_name = FALSE),
                col = cols, #c("white", "lightblue", "darkblue"), 
                # show_row_names = FALSE, 
                # show_column_names = FALSE, 
                rect_gp = gpar(col = "grey", lwd = 0.5),
                border_gp = gpar(col = "black", lty = 1), 
                ...
  ))
  
  types = sort_key$types
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(types), function(x) which(types == x)[1])-1, length(types))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # Vertical
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      })
      
      # Horizontal
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
}

SAMapAlignmentHeatmap = function(MappingTable, cols = c("white", "#fbd9d3", "#ffb09c", "red", "darkred"), add.breaks = TRUE, show_annotation_legend = FALSE, ...){
  
  # Cell type ordering
  table_order = as.vector(sapply(names(pr_palette), function(type) sapply(names(species_palette2), function(species) paste0(species, "_", type))))
  cnames = factor(colnames(MappingTable), levels = table_order)
  MappingTable = MappingTable[order(cnames), order(cnames)]
  types = str_split_fixed(colnames(MappingTable), "_", 2)[,2]
  species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
  
  # Set same species alignment scores to NA since these are never computed
  MappingTable[outer(species, species, FUN = function(x,y) x == y)] = NA
  
  # Type ordering
  row_km = factor(factor(types, levels = names(pr_palette)))
  
  rotate_and_save_image <- function(image_path, angle, output_dir) {
    image <- image_read(image_path)
    image <- image_rotate(image, angle)
    output_path <- file.path(output_dir, paste0("rotated_", basename(image_path)))
    image_write(image, output_path)
    return(output_path)
  }
  
  rotated_image_paths <- sapply(image_paths, rotate_and_save_image, angle = -90, output_dir = "../../figures/animals/black/")
  
  print(Heatmap(MappingTable, name = "Alignment\nscore", cluster_rows = FALSE, cluster_columns = FALSE, 
                width = unit(4, "in"), height = unit(4, "in"), 
                top_annotation = HeatmapAnnotation(species = anno_image(rotated_image_paths[match(str_split_fixed(colnames(MappingTable), "_", 2)[,1], names(image.files2))], 
                                                                        border = FALSE, 
                                                                        space = unit(0.1, "mm")),
                                                   `cell type` = row_km,
                                                   col = list(`cell type` = pr_palette), 
                                                   border = TRUE, 
                                                   show_annotation_name = TRUE, 
                                                   show_legend = FALSE), 
                left_annotation = rowAnnotation(species = anno_image(image.files[match(str_split_fixed(colnames(MappingTable), "_", 2)[,1], names(image.files2))], 
                                                                     border = FALSE, 
                                                                     space = unit(0.1, "mm")),
                                                `cell type` = row_km,
                                                col = list(`cell type` = pr_palette), 
                                                border = TRUE, 
                                                show_legend = show_annotation_legend,
                                                show_annotation_name = FALSE),
                col = cols, #c("white", "lightblue", "darkblue"), 
                show_row_names = FALSE, show_column_names = FALSE, 
                border_gp = gpar(col = "black", lty = 1), 
                ...
  ))
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(types), function(x) which(types == x)[1])-1, length(types))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # Vertical
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      })
      
      # Horizontal
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
}

RMSE = function(x, y){
  sqrt(mean((x-y)^2))
}

# Inspired by https://github.com/satijalab/seurat/issues/2520
PrettyUmap2 = function(object, group.by = "seurat_clusters", 
                       remove.space.x = TRUE, remove.space.y = TRUE, 
                       nbreaks = 20, geom.label = geom_text_repel, angle = 0, pad = 1.5, 
                       cols = NULL, label = TRUE, show.legend = FALSE, shuffle = TRUE, 
                       density = FALSE, alpha = 1, pt.size = 1, raster.dpi = 500,
                       centroid_fun = median, legend.point.size = 4, legend.font.size = 12,
                       ...){
  
  if(inherits(object, 'Seurat')){
    data <- as.data.frame(Embeddings(object = object[["umap"]])[(colnames(object)), c(1, 2)])
    
    # Add colors
    data$color.by <- object@meta.data[[group.by]]
  } else {
    data = object
    data$color.by = data[[group.by]]
  }
  
  # Remove empty space by binning and removing chunks with no counts
  if(remove.space.x){
    # largest x gap
    hist_x = hist(data$UMAP_1, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_x$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_index = zero_indices[which.max(counts)]
      firstZero = sum(rle$lengths[1:(longest_zero_index-1)])+1
      secondZero = firstZero + counts[which.max(counts)] - 1
      # secondZero = sum(rle$lengths[1:(longest_zero_index-2 + counts[which.max(counts)])])
      # Check 
      stopifnot(hist_x$counts[firstZero] == 0)
      stopifnot(hist_x$counts[secondZero] == 0)
      x_break1 = hist_x$mids[firstZero]
      x_break2 = hist_x$mids[secondZero]
      data[data$UMAP_1 < x_break1,"UMAP_1"] = data[data$UMAP_1 < x_break1,"UMAP_1"] + (x_break2 - x_break1)
    }
  }
  
  if(remove.space.y){
    # largest y gap
    hist_y = hist(data$UMAP_2, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_y$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_index = zero_indices[which.max(counts)]
      firstZero = sum(rle$lengths[1:(longest_zero_index-1)])+1
      secondZero = firstZero + counts[which.max(counts)] - 1
      # Check 
      stopifnot(hist_y$counts[firstZero] == 0)
      stopifnot(hist_y$counts[secondZero] == 0)
      y_break1 = hist_y$mids[firstZero]
      y_break2 = hist_y$mids[secondZero]
      data[data$UMAP_2 < y_break1,"UMAP_2"] = data[data$UMAP_2 < y_break1,"UMAP_2"] + (y_break2 - y_break1)
    }
  }
  
  # Rotate by angle
  if(angle != 0){
    data[,c("UMAP_1", "UMAP_2")] = rotate_coords(data[,c('UMAP_1', 'UMAP_2')], angle, radians = TRUE) #spdep::Rotation(data[,c("UMAP_1", "UMAP_2")], angle) %>% 
      # as.data.frame %>% 
      # setNames(c("UMAP_1", "UMAP_2"))
  }
  
  # Compute centroids for labeling
  raw_centroids = data.frame(x = sapply(unique(data$color.by), function(x) centroid_fun(data$UMAP_1[data$color.by == x])),
                             y = sapply(unique(data$color.by), function(x) centroid_fun(data$UMAP_2[data$color.by == x])),
                             NAME = unique(data$color.by))
  
  centroids <- do.call(rbind, lapply(unique(data$color.by), FUN = function(x) {
    group.data <- data[data$color.by == x, c('UMAP_1', 'UMAP_2')]
    nearest.point <- RANN::nn2(data = group.data[, 1:2], query = as.matrix(raw_centroids[raw_centroids$NAME == x, 1:2]), k = 1)$nn.idx
    data.frame(x = group.data[nearest.point, 1], 
               y = group.data[nearest.point, 2], 
               NAME = x)
  }))
  
  if(shuffle){
    data = data[sample(1:nrow(data), nrow(data)),]
  } else {
    data = data %>% arrange(desc(color.by))
  }
  
  # Plot
  plt = theme_umap(
    ggplot(data, aes(x=UMAP_1, y=UMAP_2))+
      ggrastr::rasterise(geom_point(show.legend = show.legend, aes(color = color.by), alpha = alpha, size = pt.size), dpi = raster.dpi)+
      {if(density) geom_density2d(show.legend = FALSE, alpha = 1, color = 'black', linewidth = 0.1, linetype = 'dashed', bins = 5)}+
      # {if(remove.space.x) scale_x_break(c(x_break1, x_break2))}+
      # {if(remove.space.y) scale_y_break(c(y_break1, y_break2))}+
      {if(label) geom.label(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      # {if(repel) geom_label_repel(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      {if(!is.null(cols)) scale_color_manual(values = cols)}+
      # guides(colour = guide_legend(override.aes = list(alpha = 1)))+
      scale_y_continuous(limits = c(min(data$UMAP_2)-pad, max(data$UMAP_2)+pad))+
      scale_x_continuous(limits = c(min(data$UMAP_1)-pad, max(data$UMAP_1)+pad))+
      theme(legend.text = element_text(size = legend.font.size))+
      guides(color = guide_legend(override.aes = list(size = legend.point.size, alpha = 1)))
  )
  
  plt + theme(legend.title=element_blank())
}

# Inspired by https://github.com/satijalab/seurat/issues/2520
PrettyUmap = function(object, group.by = "seurat_clusters", 
                      remove.space.x = TRUE, remove.space.y = TRUE, 
                      nbreaks = 20, geom.label = geom_text_repel, angle = 0, pad = 0.5, 
                      cols = NULL, label = TRUE, show.legend = FALSE, shuffle = TRUE, 
                      density = FALSE, alpha = 0.03, pt.size = 3, raster.dpi = 500,
                      ...){
  
  if(inherits(object, 'Seurat')){
    data <- as.data.frame(Embeddings(object = object[["umap"]])[(colnames(object)), c(1, 2)])
    
    # Add colors
    data$color.by <- object@meta.data[[group.by]]
  } else {
    data = object
    data$color.by = data[[group.by]]
  }
  
  # Remove empty space by binning and removing chunks with no counts
  if(remove.space.x){
    # largest x gap
    hist_x = hist(data$UMAP_1, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_x$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_stretch = 
        firstZero = sum(rle$lengths[1:(which.max(rle$lengths)-1)])+1
      secondZero = sum(rle$lengths[1:which.max(rle$lengths)])
      x_break1 = hist_x$mids[firstZero]
      x_break2 = hist_x$mids[secondZero]
      data[data$UMAP_1 < x_break1,"UMAP_1"] = data[data$UMAP_1 < x_break1,"UMAP_1"] + (x_break2 - x_break1)
    }
  }
  
  if(remove.space.y){
    # largest y gap
    hist_y = hist(data$UMAP_2, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_y$counts)
    firstZero = sum(rle$lengths[1:(which.max(rle$lengths)-1)])
    secondZero = sum(rle$lengths[1:which.max(rle$lengths)])
    y_break1 = hist_y$mids[firstZero]
    y_break2 = hist_y$mids[secondZero]
    data[data$UMAP_2 < y_break1,"UMAP_2"] = data[data$UMAP_2 < y_break1,"UMAP_2"] + (y_break2 - y_break1)
  }
  
  # Rotate by angle
  data[,c("UMAP_1", "UMAP_2")] = spdep::Rotation(data[,c("UMAP_1", "UMAP_2")], angle) %>% 
    as.data.frame %>% 
    setNames(c("UMAP_1", "UMAP_2"))
  
  # Compute centroids for labeling
  centroids = data.frame(x = sapply(unique(data$color.by), function(x) mean(data$UMAP_1[data$color.by == x])),
                         y = sapply(unique(data$color.by), function(x) mean(data$UMAP_2[data$color.by == x])), 
                         NAME = unique(data$color.by))
  
  if(shuffle){
    data = data[sample(1:nrow(data), nrow(data)),]
  } else {
    data = data %>% arrange(desc(color.by))
  }
  
  # Plot
  plt = theme_umap(
    ggplot(data, aes(x=UMAP_1, y=UMAP_2))+
      ggrastr::rasterise(geom_point(show.legend = show.legend, aes(color = color.by), alpha = alpha, size = pt.size), dpi = raster.dpi)+
      {if(density) geom_density2d(show.legend = FALSE, alpha = 1, color = 'black', linewidth = 0.1, linetype = 'dashed', bins = 5)}+
      # {if(remove.space.x) scale_x_break(c(x_break1, x_break2))}+
      # {if(remove.space.y) scale_y_break(c(y_break1, y_break2))}+
      {if(label) geom.label(data = centroids, aes(label = NAME, x = x, y = y), point.size = NA, ...)}+
      # {if(repel) geom_label_repel(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      {if(!is.null(cols)) scale_color_manual(values = cols)}+
      guides(colour = guide_legend (override.aes = list(alpha = 1)))+
      scale_y_continuous(limits = c(min(data$UMAP_2)-pad, max(data$UMAP_2)+pad))+
      scale_x_continuous(limits = c(min(data$UMAP_1)-pad, max(data$UMAP_1)+pad))
  ) 
  
  plt + theme(legend.title=element_blank())
}

DrawSankey4 = function(filepath, min.value = 0.1, species = c('ze', 'ch', 'li', 'op', 'rn', 'hs'), axis.width = 0.1){
  library(ggforce)
  
  species1 = 'ze'
  species2 = 'ch'
  species3 = 'li'
  species4 = 'op'
  species5 = 'rn'
  species6 = 'hs'
  
  color.mapping = pr_palette
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  # print(head(MappingTable))
  MappingTable = rbind(MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species2) & startsWith(as.character(MappingTable$Var2), species3),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species3) & startsWith(as.character(MappingTable$Var2), species4),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species4) & startsWith(as.character(MappingTable$Var2), species5),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species5) & startsWith(as.character(MappingTable$Var2), species6),]
  )
  
  # print(head(MappingTable))
  # MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
  #                                         factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Ordering of strata
  # MappingTable$Var1 = factor(MappingTable$Var1, levels = paste0(species1, '_', names(pr_palette)))
  # MappingTable$Var2 = factor(MappingTable$Var2, levels = paste0(species2, '_', names(pr_palette)))
  # MappingTable$Var3 = factor(MappingTable$Var3, levels = paste0(species3, '_', names(pr_palette)))
  # MappingTable$Var4 = factor(MappingTable$Var4, levels = paste0(species4, '_', names(pr_palette)))
  # MappingTable$Var5 = factor(MappingTable$Var5, levels = paste0(species5, '_', names(pr_palette)))
  # MappingTable$Var6 = factor(MappingTable$Var6, levels = paste0(species6, '_', names(pr_palette)))
  MappingTable = gather_set_data(MappingTable, 1:2)
  
  print((MappingTable))
  MappingTable$names = factor(str_split_fixed(MappingTable$y, '_', 2)[,2], levels = names(pr_palette))
  
  # Color nodes
  colors <- rep(pr_palette, 6)
  names(colors) = c(paste0(species1, '_', names(pr_palette)), 
                    paste0(species2, '_', names(pr_palette)), 
                    paste0(species3, '_', names(pr_palette)), 
                    paste0(species4, '_', names(pr_palette)), 
                    paste0(species5, '_', names(pr_palette)), 
                    paste0(species6, '_', names(pr_palette)))
  print(colors)
  MappingTable = MappingTable[MappingTable$value > min.value,]
  # return(MappingTable)
  
  # Plot
  sn = ggplot(MappingTable, aes(x, id = id, split = names, value = value)) +
    geom_parallel_sets(axis.width = axis.width, fill = 'grey', alpha = 0.4) +
    geom_parallel_sets_axes(aes(fill = y), axis.width = axis.width, color = "black") +
    geom_parallel_sets_labels(colour = 'black',
                              angle = 0
                              # hjust = c(rep(1, length(unique(MappingTable$Var1))), rep(0, length(unique(MappingTable$Var2)))),
                              # nudge_x = c(rep(-0.08, length(unique(MappingTable$Var1))), rep(0.08, length(unique(MappingTable$Var2))))
    ) +
    scale_fill_manual(values = colors) +
    theme_void() +
    NoLegend() +
    scale_x_discrete(expand = expansion(add = 0.5))
  
  return(sn)
}

DrawSankey3 = function(input, min.value = 0.1, species1 = 'ch', species2 = 'li', axis.width = 0.07, label = TRUE, text_nudge = 0.07){
  library(ggforce)
  
  color.mapping = pr_palette
  MappingTable = reshape2::melt(input) #as.matrix(read.csv(filepath, row.names = 1)))
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                          factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Ordering of strata
  MappingTable$Var1 = factor(MappingTable$Var1, levels = paste0(species1, '_', names(pr_palette)))
  MappingTable$Var2 = factor(MappingTable$Var2, levels = paste0(species2, '_', names(pr_palette)))
  MappingTable = gather_set_data(MappingTable, 1:2)
  
  # MappingTable$Var1 = str_split_fixed(MappingTable$Var1, '_', 2)[,2]
  # MappingTable$Var2 = str_split_fixed(MappingTable$Var2, '_', 2)[,2]
  MappingTable$names = factor(str_split_fixed(MappingTable$y, '_', 2)[,2], levels = names(pr_palette))
  
  # Color nodes
  colors <- rep(pr_palette, 2)
  names(colors) = c(paste0(species1, '_', names(pr_palette)), 
                    paste0(species2, '_', names(pr_palette)))
  
  MappingTable = MappingTable[MappingTable$value > min.value,]
  
  # Plot
  sn = ggplot(MappingTable, aes(x, id = id, split = names, value = value)) +
    geom_parallel_sets(axis.width = axis.width, fill = 'grey', alpha = 0.4) +
    geom_parallel_sets_axes(aes(fill = y), axis.width = axis.width, color = "black") +
    {if(label) geom_parallel_sets_labels(colour = 'black', 
                                         angle = 0,
                                         hjust = c(rep(1, length(unique(MappingTable$Var1))), rep(0, length(unique(MappingTable$Var2)))),
                                         nudge_x = c(rep(-text_nudge, length(unique(MappingTable$Var1))), rep(text_nudge, length(unique(MappingTable$Var2)))))} +
    scale_fill_manual(values = colors) +
    theme_void() +
    NoLegend() +
    scale_x_discrete(expand = expansion(add = 0.5))
  
  return(sn)
}

GetTypeList = function(gene_pairs, PHOTORECEPTORS, species_ids){
  # gene.summary = list()
  type.list = list()
  
  # This recovers all 335 entries
  for(i in 1:(length(species_ids)-1)){
    type.list[[ species_ids[i] ]] = list()
    for(k in 1:(length(PHOTORECEPTORS))){
      type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]] = list()
      for(j in (i+1):length(species_ids)){
        type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]][[ species_ids[j] ]] = list()
        for(l in (1):length(PHOTORECEPTORS)){
          # element = 1
          column.name = paste0(species_ids[i], "_", PHOTORECEPTORS[k], ".", species_ids[j], "_", PHOTORECEPTORS[l])
          if(column.name %in% names(gene_pairs)){
            genes.1 = ExtractString(gene_pairs[[column.name]], paste0(species_ids[i], "_"), ";")
            genes.1 = genes.1[genes.1 != ""]
            genes.2 = ExtractString(gene_pairs[[column.name]], paste0(species_ids[j], "_"))
            genes.2 = genes.2[genes.2 != ""]
            type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]][[ species_ids[j] ]][[ PHOTORECEPTORS[l] ]] = setNames(data.frame(genes.1, genes.2), names(species_ids)[c(i,j)])
            # element = element + 1
          }
        }
      }
    }
  }
  
  return(type.list)
}

DrawSankey2 = function(filepath, min.value = 0.1, species1, species2){
  
  color.mapping = pr_palette
  # color.mapping = c("grey", "violet", "blue", "green", "red", "magenta", "cyan", 'black', 'brown') %>% setNames(c('rod', 'UV', 'blue', 'green', 'red', 'principle', 'accessory', 'AC', 'HC'))
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                          factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  alignment = list()
  alignment$nodes = data.frame(name = unique(c(MappingTable$Var1, MappingTable$Var2)))
  nodeIDs = alignment$nodes$name
  names(nodeIDs) = as.numeric(rownames(alignment$nodes))-1 # Zero index
  
  # Filtration
  alignment$links = MappingTable
  
  # Convert to index
  alignment$links$Var1 = as.numeric(names(nodeIDs[match(alignment$links$Var1, nodeIDs)]))
  alignment$links$Var2 = as.numeric(names(nodeIDs[match(alignment$links$Var2, nodeIDs)]))
  alignment$links = alignment$links[!is.na(alignment$links$value),]
  alignment$links = alignment$links[alignment$links$value >= min.value,]
  # alignment$links$value[is.na(alignment$links$value)] = 0
  rownames(alignment$links) = 1:nrow(alignment$links)
  
  # Color nodes
  alignment$nodes$node_color = str_split_fixed(alignment$nodes$name, "_", 2)[,2]
  colors <- paste(unique(color.mapping[match(alignment$nodes$node_color, names(color.mapping))]), collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  alignment$nodes$display = ""
  
  # Plot
  sn = sankeyNetwork(Links = alignment$links, 
                     Nodes = alignment$nodes, Source = 'Var1',
                     Target = 'Var2', Value = 'value', NodeID = 'node_color',
                     fontSize = 12, nodeWidth = 30, NodeGroup = "node_color", 
                     colourScale = colorJS, fontFamily = "arial")
  
  saveNetwork(sn, "sn.html")
  webshot::webshot("sn.html", gsub(".csv", ".png", filepath), vwidth = 300, vheight = 400, zoom = 2)
  
  return(sn)
}

DrawSankey = function(filepath, min.value = 0.1, species1, species2){
  
  color.mapping = c("purple", "blue", "green", "red", "grey", "cyan", "magenta", 'black', 'brown') %>% setNames(c('UV', 'blue', 'green', 'red', 'rod', 'accessory', 'principle', 'AC', 'HC'))
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  alignment = list()
  alignment$nodes = data.frame(name = unique(MappingTable$Var2))
  nodeIDs = alignment$nodes$name
  names(nodeIDs) = as.numeric(rownames(alignment$nodes))-1 # Zero index
  
  # Filtration
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  alignment$links = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                             factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Convert to index
  alignment$links$Var1 = as.numeric(names(nodeIDs[match(alignment$links$Var1, nodeIDs)]))
  alignment$links$Var2 = as.numeric(names(nodeIDs[match(alignment$links$Var2, nodeIDs)]))
  alignment$links = alignment$links[!is.na(alignment$links$value),]
  alignment$links = alignment$links[alignment$links$value >= min.value,]
  # alignment$links$value[is.na(alignment$links$value)] = 0
  rownames(alignment$links) = 1:nrow(alignment$links)
  
  # Color nodes
  alignment$nodes$node_color = str_split_fixed(alignment$nodes$name, "_", 2)[,2]
  colors <- paste(unique(color.mapping[match(alignment$nodes$node_color, names(color.mapping))]), collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  alignment$nodes$display = ""
  
  # Plot
  sn = sankeyNetwork(Links = alignment$links, 
                     Nodes = alignment$nodes, Source = 'Var1',
                     Target = 'Var2', Value = 'value', NodeID = 'node_color',
                     fontSize = 12, nodeWidth = 30, NodeGroup = "node_color", 
                     colourScale = colorJS, fontFamily = "arial")
  
  saveNetwork(sn, "sn.html")
  webshot::webshot("sn.html", gsub(".csv", ".png", filepath), vwidth = 300, vheight = 400, zoom = 2)
  
  return(sn)
}

RenameFeatures = function(object, new.names){
  
  if(length(unique(new.names)) != length(new.names)) message("Found ", length(new.names)-length(unique(new.names))," duplicates...making unique!")
  
  rownames(object@assays$RNA@counts) <- make.unique((new.names))
  rownames(object@assays$RNA@data) <- make.unique((new.names))
  rownames(object@assays$RNA@meta.features) <- make.unique((new.names))
  
  # Would need to match order for this
  
  # if('integrated' %in% names(objectList$Lizard@assays)){
  #   rownames(object@assays$integrated@counts) <- make.unique((new.names))
  #   rownames(object@assays$integrated@data) <- make.unique((new.names))
  #   rownames(object@assays$integrated@meta.features) <- make.unique((new.names))
  # }
    
  # rownames(object@assays$RNA@scale.data)<- make.unique((new.names))
  
  return(object)
}

# VlnPlot2 = function(object){
#   VlnPlot(object, ...) + 
#     NoLegend() + 
#     theme(axis.title)
# }

SelectFeatures = function(object, nfeatures = 1000, min.pct.expressed = 10, group.by = "seurat_clusters", assay = "RNA", return.pct = FALSE){
  
  object = FindVariableFeatures(object, selection.method = "vst", nfeatures = nrow(object@assays[[assay]]@data), verbose = FALSE)
  
  # Filter out lowly expressed genes
  if(min.pct.expressed > 0){
    pct.expressed = PercentageExpressed(object, features = object@assays[[assay]]@var.features, group.by = group.by)
    pct.expressed$max.exp = apply(pct.expressed, 1, max)
    if(return.pct) return(pct.expressed)
    features.use = object@assays[[assay]]@var.features[pct.expressed$max.exp >= min.pct.expressed]
  } else {
    features.use = object@assays[[assay]]@var.features
  }
  
  features.use = head(features.use, nfeatures)
  message("Using top ", length(features.use)," variable features expressed in at least ", min.pct.expressed, "% cells of at least one cluster")
  
  return(features.use)
}

PseudoBulkPCA = function(object, assay = "RNA", group.by = "seurat_clusters", features = NULL, title = NULL, 
                         return.plot = TRUE, min.pct.expressed = 0, nfeatures = 2000, axes = c(1,2), scale = FALSE, binary = NULL, 
                         size = 1, shape = 1, cols = NULL, label = TRUE, annotate.by = NULL){
  
  if(inherits(object, "Seurat")){
    
    if(!is.null(features)){
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features)[[assay]]))
    } else if(min.pct.expressed == 0){
      object = FindVariableFeatures(object, nfeatures = nfeatures)
      features.use = object@assays[[assay]]@var.features
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features.use)[[assay]]))
    } else {
      # Feature selection to avoid noisy lowly expressed genes that get amplified after scaling
      features.use = SelectFeatures(object, nfeatures = nfeatures, min.pct.expressed = min.pct.expressed, group.by = group.by)
      
      # Get average expression
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features.use)[[assay]]))
    }
  } else {
    avg_exp = object
  }
  
  print(dim(avg_exp))
  
  if(!is.null(binary)){
    avg_exp_old = avg_exp
    avg_exp[avg_exp_old >= binary] = 1
    avg_exp[avg_exp_old < binary] = 0
  }
  
  # Compute principle components
  res.pca <- prcomp(t(avg_exp), scale = scale)
  if(!return.plot) return(res.pca)
  
  var_exp = (res.pca$sdev^2)/sum(res.pca$sdev^2)*100
  
  # Plot data projected onto PC1 and PC2
  df = as.data.frame(res.pca$x)
  df$label = rownames(df)
  if(!is.null(annotate.by)) {
    df$group = Metadata(object, group.by, annotate.by)[,annotate.by]
  } else {
    df$group = df$label
  }
  
  TitlePlot(
    ggplot(df, aes(x = PC1, y = PC2, label = label, color = group))+
      geom_point(size = size, shape = shape)+
      {if(label) geom_text_repel(show.legend = FALSE)}+
      {if(!is.null(cols)) scale_color_manual(values = cols)}+
      labs(x = paste0("PC1 (", signif(var_exp[1], 2), "%)"), y = paste0("PC2 (", signif(var_exp[2], 2), "%)"))+
      theme_dario()+
      theme(legend.key=element_rect(fill="white")),
    
    # factoextra::fviz_pca_ind(res.pca,
    #                          axes = axes,
    #                          # habillage=
    #                          # col.ind = "cos2", # Color by the quality of representation
    #                          # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    #                          repel = TRUE,     # Avoid text overlapping
    #                          axes.linetype=NA
    # ) + theme_dario(),
    title = title
  )
  
}

# SAMap genes heatmap!
SAMapHeatmap2 = function(objectList, 
                        gene.key, 
                        type_palette, 
                        species_palette,
                        species.use = c("Chicken", "Zebrafish"), 
                        types.use = c("red", "green", "blue", "UV", "rod"), 
                        min.z.score = -2,
                        max.z.score = 2, 
                        rotate = FALSE,
                        dotplot = FALSE, 
                        dot.scale.factor = 0.05,
                        col_fun = NULL,
                        names.show = NULL,
                        ...){
  
  args = list(...)
  
  # Set col_fun
  if(is.null(col_fun)){
    col_fun = circlize::colorRamp2(c(min.z.score, 0, max.z.score), c("blue", "white", "red"))
  }
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    
    norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    colnames(scaled.expr) = paste0(species, ' ', colnames(scaled.expr))
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < min.z.score] = min.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Percentage expression
  pct.expr = lapply(species.use, function(species) {
    pct.expr = PercentageExpressed(subset(objectList[[species]], annotated %in% types.use), features = setdiff(unique(gene.key[[species]]), ''), group.by = 'annotated')
    pct.expr = pct.expr[match(gene.key[[species]], rownames(pct.expr)),]
    colnames(pct.expr) = paste0(species, ' ', colnames(pct.expr))
    pct.expr
  })
  
  # Combine 
  if(is.null(names.show)) {
    gene.names = apply(gene.key, 1, function(x) paste0(x, collapse = "/"))
  } else {
    gene.names = make.unique(gene.key[[names.show]])
  }
  
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  combined.pct = do.call(cbind, pct.expr)
  rownames(combined.pct) = gene.names
  
  # Ordering matrix
  types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(y, " ", x)), 1, function(x) x))
  combined.expr = combined.expr[, order(factor(colnames(combined.expr), levels = types.order))]
  combined.pct = combined.pct[, order(factor(colnames(combined.pct), levels = types.order))]
  
  # Ordering for annotations
  types = str_split_fixed(colnames(combined.expr), " ", 2)[,2]
  species = str_split_fixed(colnames(combined.expr), " ", 2)[,1]
  row_split = paste0(gene.key$row_annotation, ' markers')
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(rotate){
    # Transpose
    combined.expr = t(combined.expr)
    combined.pct = t(combined.pct)
  } 
  
  if(dotplot){
    min.z.score = -1
    dot_col = function(species, expression) circlize::colorRamp2(c(min.z.score, max.z.score), c('white', species_palette[[species]]))(expression)
    cell_fun = function(j, i, x, y, w, h, fill){
      grid.circle(x=x,
                  y=y,
                  r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                  gp=gpar(fill = dot_col(species[[j]], combined.expr[i, j]), col = NA))}
                  # gp=gpar(fill = col_fun(combined.expr[i, j]), col = NA))}
    rect_gp = gpar(type = "none")
    
    # Set this for color legend
    col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c('white', 'grey25'))
    
    # Create a dot legend
    lgd = Legend(labels = seq(0,100, by = 25), title = "Percentage\nexpressed",
                 graphics = list(
                   function(x, y, w, h) grid.circle(x, y, r=unit(0/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(25/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(50/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(75/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(100/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
                 ))
  } else {
    cell_fun = NULL
    rect_gp = gpar(col = NA)
    lgd = NULL
  }
  
  # Plot
  if(rotate){
    draw(Heatmap(combined.expr,
                  # name = "mat",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  heatmap_legend_param = list(title = "Scaled\nexpression"),
                  column_split = factor(row_split, levels = unique(row_split)), 
                  row_split = factor(types, levels = unique(types)),
                  row_title = NULL,
                  border = TRUE, 
                  cell_fun = cell_fun,
                  rect_gp = rect_gp,
                  left_annotation = rowAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))]),
                                                                           labels = unique(types), 
                                                                           labels_gp = gpar(col = "black", fontsize = 10)),
                                                  Species = species,
                                                  # annotation_name_side = "bottom", 
                                                  show_annotation_name = FALSE,
                                                  show_legend = args$show_heatmap_legend, 
                                                  border = TRUE, 
                                                  col = list(`Cell type` = type_palette, Species = species_palette
                                                  )), 
                  ...), annotation_legend_list = lgd)
  } else {
    draw(Heatmap(combined.expr,
                  # name = "mat",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  row_split = factor(row_split, levels = unique(row_split)),
                  column_split = factor(types, levels = unique(types)),
                  column_title = NULL,
                  border = TRUE, 
                  heatmap_legend_param = list(title = "Scaled\nexpression"),
                  cell_fun = cell_fun,
                 # cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", combined.pct[i, j]), x, y, gp = gpar(fontsize = 10))},
                  rect_gp = rect_gp,
                  top_annotation = HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))]),
                                                                              labels = unique(types), 
                                                                              labels_gp = gpar(col = "black", fontsize = 10)),
                                                     Species = species,
                                                     show_annotation_name = FALSE,
                                                     # annotation_name_side = 'none', 
                                                     show_legend = args$show_heatmap_legend, 
                                                     border = TRUE, 
                                                     col = list(`Cell type` = type_palette, Species = species_palette
                                                     )), 
                  ...), annotation_legend_list = lgd)
  }
}

# SAMap genes heatmap!
SAMapHeatmap = function(objectList, 
                        gene.key, 
                        species.use = c("Chicken", "Zebrafish"), 
                        types.use = c("red", "green", "blue", "UV", "rod"), 
                        max.z.score = 2, 
                        rotate = FALSE,
                        ...){
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data")$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    colnames(scaled.expr) = paste0(colnames(scaled.expr), "_", species)
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < -max.z.score] = -max.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Combine 
  gene.names = do.call(cbind, lapply(scaled.expr, function(x) rownames(x)))
  gene.names = apply(gene.names, 1, function(x) paste0(x, collapse = "/"))
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  # Ordering matrix
  types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(x, "_", y)), 1, function(x) x))
  combined.expr = combined.expr[, order(factor(colnames(combined.expr), levels = types.order))]
  
  # Ordering for annotations
  types = (str_split_fixed(colnames(combined.expr), "_", 2)[,1])
  species = str_split_fixed(colnames(combined.expr), "_", 2)[,2]
  row_split = gene.key$row_annotation
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(rotate){
    print(Heatmap(t(combined.expr),
                  name = "mat",
                  col = colorRamp2(c(-max.z.score, 0, max.z.score), c("blue", "white", "red")),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  heatmap_legend_param = list(title = ""),
                  column_split = factor(row_split, levels = unique(row_split)), 
                  border = TRUE, 
                  left_annotation = rowAnnotation(species = anno_image(image.files[match(species, names(image.files))], 
                                                                       border = FALSE, 
                                                                       space = unit(0, "mm")),
                                                  `cell type` = str_split_fixed(colnames(combined.expr), "_", 2)[,1],
                                                  annotation_name_side = "bottom", 
                                                  show_legend = FALSE, 
                                                  border = TRUE, 
                                                  col = list(`cell type` = type_palette
                                                  )), 
                  ...))
    
    # Add horizontal bars separating rows
    types = rev(types)
    breaks = (sapply(unique(types), function(x) which(types == x)[1])-1)[-1]
    for(index in breaks) {
      for(slice in seq_along(unique(row_split))){
        x_coord = index/length(types)
        
        # Horizontal
        decorate_heatmap_body("mat", column_slice = slice, {
          grid.lines(c(0, 1), c(x_coord, x_coord), gp = gpar(lty = 1, lwd = 1))
        })
      }
    }
  } else {
    print(Heatmap(combined.expr,
                  name = "mat",
                  col = colorRamp2(c(-max.z.score, 0, max.z.score), c("blue", "white", "red")),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  row_split = factor(row_split, levels = unique(row_split)), 
                  border = TRUE, 
                  heatmap_legend_param = list(title = ""),
                  # rect_gp = gpar(col = "black", lwd = 1),
                  # left_annotation = rowAnnotation(#row = anno_textbox(gene.key$row_annotation, gene.key$row_annotation), 
                  #                     row = anno_block(labels = unique(gene.key$row_annotation), 
                  #                                labels_gp = gpar(col = "white", fontsize = 10))),
                  top_annotation = HeatmapAnnotation(species = anno_image(image.files[match(species, names(image.files))], 
                                                                          border = FALSE, 
                                                                          space = unit(0, "mm")),
                                                     `cell type` = str_split_fixed(colnames(combined.expr), "_", 2)[,1],
                                                     annotation_name_side = "left", 
                                                     show_legend = FALSE, 
                                                     border = TRUE, 
                                                     # gp = gpar(col = "black"), 
                                                     col = list(`cell type` = pr_palette
                                                     )), 
                  ...))
    
    # Add horizontal bars separating rows
    breaks = (sapply(unique(types), function(x) which(types == x)[1])-1)[-1]
    for(index in breaks) {
      for(slice in seq_along(unique(row_split))){
        x_coord = index/length(types)
        
        # Vertical
        decorate_heatmap_body("mat", row_slice = slice, {
          grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
        })
      }
    }
  }
}

# From Seurat DotPlot function
DarPlot = function(data, dot.scale = 6, scale.min = NA, scale.max = NA){
  plot <- ggplot(data = data, mapping = aes_string(x = "features.plot", y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
    scale_color_gradient(low = "lightgrey", high = "blue") + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme_cowplot() + coord_flip() + RotatedAxis() + theme(axis.title = element_blank())
  
  plot
}

TopNDEGs = function(object, group.by = "annotated", n = 10, avg_log2FC_cutoff = 0.25, p_val_adj_cutoff = 0.05){
  
  if(inherits(object, 'Seurat')){
    # de_table = wilcoxauc(object, group_by = group.by) %>% 
    #   setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
    de_table = FindAllMarkersFast(object, group.by = group.by, avg_log2FC_cutoff = avg_log2FC_cutoff, p_val_adj_cutoff = p_val_adj_cutoff)
  } else {
    de_table = object
  }
  
  de_table %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff) %>%
    arrange(-avg_log2FC) %>% 
    slice_head(n = n) %>%
    ungroup() -> top
  
  if(is.factor(object@meta.data[[group.by]])){
    # top = top %>% arrange(factor(rlang::sym(c("annotated")), levels = levels(object@meta.data[[group.by]])))
    top$cluster = factor(top$cluster, levels = levels(object@meta.data[[group.by]]))
    top = top[order(top$cluster),]
  }
  
  top %>% as.data.frame 
}

CelltypeProportionBarplot = function(object, x = "annotated", y = "animal", return.table = FALSE, show.all = FALSE, normalize = TRUE){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  # Remove unused factors
  if(is.factor(data[[x]])) data[[x]] = factor(factor(data[[x]], levels = (levels(data[[x]]))))
  
  tabulation = table(data[[y]], data[[x]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  if(show.all){
    p = ggplot(melted, aes(y = Proportion, x = Type, fill = Sample))+
      geom_bar(stat = "identity", position="dodge") + 
      theme_bw() + 
      scale_fill_discrete(name = "Sample")+
      ylab("Proportion") 
  } else {
    p = ggbarplot(melted, y = "Proportion", x = "Type", fill = "Type", add = c("mean_se", "jitter")) + 
      NoLegend() +
      scale_y_continuous(expand = expansion(mult = c(0, .1)))
    # rremove("xlab")
  }
  if(return.table) {
    if(normalize) {
      normalized.tabulation
    } else {
      tabulation
    }
  } else {
    p
  }
}

PercentageExpressed = function(object, features, group.by = NULL){
  
  if(!is.null(group.by)){
    obj.list = SplitObject(object, split.by = group.by)
    pct.exp = as.data.frame(do.call(cbind, lapply(obj.list, function(obj){
      datExpr = obj@assays$RNA@counts[features, , drop = FALSE]
      pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
    })))
    
    # Sort by levels of group.by
    if(is.factor(object@meta.data[[group.by]])) pct.exp = pct.exp[,levels(object@meta.data[[group.by]])]
  } else {
    datExpr = object@assays$RNA@counts[features, , drop = FALSE]
    pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
  }
  
  return(pct.exp)
  
  # tabulation = (object@assays$RNA@counts[feature,] > 0) %>% table
  # pct.exp = as.numeric((tabulation[["TRUE"]]/sum(tabulation)))*100
  # return(pct.exp)
}

# From https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

BlastOrthologyTable = function(file1 = "maps/zech/ch_to_ze.txt", 
                               file2 = "maps/zech/ze_to_ch.txt", 
                               keyfile1 = "Keys/Chicken_key.csv",
                               keyfile2 = "Keys/Zebrafish_key.csv", 
                               one2one = TRUE){
  
  orthology_file1 = fread(file1) %>% as.data.frame
  orthology_file2 = fread(file2) %>% as.data.frame
  key1 = fread(keyfile1, header = T) %>% as.data.frame
  key2 = fread(keyfile2, header = T) %>% as.data.frame
  
  # Convert transcripts to genes
  orthology_file1$V1 = key1$symbol[match(orthology_file1$V1, key1[,2])]
  orthology_file1$V2 = key2$symbol[match(orthology_file1$V2, key2[,2])]
  orthology_file2$V1 = key2$symbol[match(orthology_file2$V1, key2[,2])]
  orthology_file2$V2 = key1$symbol[match(orthology_file2$V2, key1[,2])]
  
  if(!one2one) return(list(orthology_file1, orthology_file2))
  
  # Make unique 
  # orthology_file1 = orthology_file1 %>% unique
  # orthology_file2 = orthology_file2 %>% unique
  
  # Get best 1 to 2 blast
  best1 = orthology_file1 %>% 
    filter(V11 < 1e-6) %>%
    dplyr::select(V1, V2, V12) %>% 
    group_by(V1) %>%
    filter(V12 == max(V12, na.rm=TRUE)) %>% 
    unique
  
  best2 = orthology_file2 %>% 
    filter(V11 < 1e-6) %>%
    dplyr::select(V1, V2, V12) %>% 
    group_by(V1) %>%
    filter(V12 == max(V12, na.rm=TRUE)) %>% 
    unique
  
  # Binary masks
  matrix1 = reshape2::acast(best1, V1 ~ V2, value.var = "V12")
  matrix1[!is.na(matrix1)] = 1
  matrix1[is.na(matrix1)] = 0
  
  matrix2 = reshape2::acast(best2, V1 ~ V2, value.var = "V12")
  matrix2[!is.na(matrix2)] = 1
  matrix2[is.na(matrix2)] = 0
  
  # Subset to common genes
  gene1 = intersect(rownames(matrix1), colnames(matrix2))
  gene2 = intersect(rownames(matrix2), colnames(matrix1))
  final1 = matrix1[gene1, gene2]
  final2 = t(matrix2[gene2, gene1])
  
  # Find bidirectional blast hits
  A = final1*final2
  
  # Remove genes with no ortholog
  A = A[rowSums(A) > 0,]
  A = A[,colSums(A) > 0]
  
  # Underscores to dashes for compatibility with seurat
  colnames(A) = gsub("_", "-", colnames(A))
  rownames(A) = gsub("_", "-", rownames(A))
  
  # Get ambiguous mappings
  dupl_rows = apply(A, 1, function(x) sum(x) > 1)
  dupl_cols = apply(A, 2, function(x) sum(x) > 1)
  
  # Remove ambiguous orthologs
  new = A[!dupl_rows,!dupl_cols]
  
  # Remove rows/columns with zeroes
  zero_rows = apply(new, 1, function(x) sum(x) == 0)
  zero_cols = apply(new, 2, function(x) sum(x) == 0)
  new2 = new[!zero_rows, !zero_cols]
  
  return(new2)
}

DotPlot2 = function(object, coord.flip = TRUE, max.pct = 100, ...){
  DotPlot(object, col.min = -1, col.max = 2, ...) + 
    {if(coord.flip) coord_flip()}+ 
    RotatedAxis() + 
    theme(axis.title = element_blank())+
    # labs(color = 'Average\nexpression', size = 'Percent\nexpressed')+
    guides(color = guide_colorbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed"))+
    scale_color_gradient(name = 'Averageexpression', low = "lightgrey", high = "#584B9FFF", limits = c(-1,2))+
    scale_radius(name = 'Percentexpressed', limits = c(0,max.pct), range = c(0, 6))
  # scale_x_discrete(expand = c(0,0))+
  # scale_y_discrete(expand = c(0,0))+
  # theme(panel.grid = element_line(color = rgb(235, 235, 235, 100, 
  #                                             maxColorValue = 255), 
  #                                 linewidth = 1, 
  #                                 linetype = 1))
}

SeuratToH5ad2 = function(object, filepath, counts.only = TRUE){
  
  # Turn to factor
  if('annotated' %in% colnames(object@meta.data)) object$annotated = as.character(object$annotated)
  
  if(counts.only){
    new_object = CreateSeuratObject(object@assays$RNA@counts)
    new_object = TransferMetadata(object, new_object)
    object = new_object
  }
  
  tempfile = gsub("h5ad", "h5Seurat", filepath)
  SeuratDisk::SaveH5Seurat(object, filename = tempfile, overwrite = TRUE)
  SeuratDisk::Convert(tempfile, dest = "h5ad", overwrite = TRUE)
  file.remove(tempfile)
  
  return(object)
}

SeuratToH5ad = function(object, filepath, genes.remove = NULL, types.remove = NULL, downsample = 100, seed = 12345){
  
  # Unfactor as this causes issues downstream
  object$annotated = as.character(object$annotated)
  
  if(!is.null(types.remove)){
    object = subset(object, annotated %in% types.remove, invert = TRUE)
  }
  if(!is.null(downsample)) object = DownsampleSeurat(object, group.by = "annotated", size = downsample, seed = seed)
  message("Some cell ids: ", paste0(head(Cells(object)), collapse = ", "))
  annotated = object$annotated
  animal = object$animal
  object = CreateSeuratObject(object@assays$RNA@counts[!rownames(object@assays$RNA@counts) %in% genes.remove,])
  object$annotated = annotated
  object$animal = animal
  print(table(object$annotated))
  
  tempfile = gsub("h5ad", "h5Seurat", filepath)
  SeuratDisk::SaveH5Seurat(object, filename = tempfile, overwrite = TRUE)
  SeuratDisk::Convert(tempfile, dest = "h5ad", overwrite = TRUE)
  file.remove(tempfile)
  
  return(object)
}

ExtractString = function(vector, before = NULL, after = NULL){
  if(is.null(after)){
    return(str_split_fixed(vector, before, 2)[,2])
  } else if(is.null(before)){
    return(str_split_fixed(vector, after, 2)[,1])
  } else {
    return(str_split_fixed(str_split_fixed(vector, after, 2)[,1], before, 2)[,2])
  }
}

UpdateCellClass = function(object, annotation = major_annotation){
  # Convert to major cell class
  object$cell_class = as.character(object$cell_class)
  object$cell_class[object$cell_class == "GabaAC" | object$cell_class == "GlyAC"] = "AC"
  object$cell_class[object$cell_class == "BP"] = "BC"
  
  # Everything else goes to other
  object$cell_class[!object$cell_class %in% Annotation(annotation)] = "Other"
  # object$cell_class[object$cell_class == "MicroG"] = "Other"
  object$cell_class = factor(object$cell_class, unique(Annotation(annotation)))
  
  return(object)
}

ProportionBarplot = function(data, x, y, return.table = FALSE){
  tabulation = table((data[[x]]), (data[[y]]))
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation)
  p = ggplot(melted, aes(y = value, x = Var2, fill = Var1))+
    geom_bar(stat = "identity", position="dodge") + 
    rremove("xlab") + ylab("Proportion")
  if(return.table) {
    normalized.tabulation
  } else {
    p + theme_dario()
  }
}

PrettyBoxplot = function(data){
  ggplot(reshape2::melt(data), aes(x = variable, y = value, color = variable))+
    geom_boxplot(linetype = "dashed", outlier.shape = NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5) +
    # ylab(ylab)+
    # ggtitle("Mixing of species")+
    theme_cowplot()+
    theme(legend.position = "none", axis.title.x = element_blank())
}

dtColors = function(n = 435, omit = c('white', 'ivory', 'grey60', 'darkgrey', 'lightyellow', 'lightcyan', 'floralwhite', 'lightcyan1')){
  colors = standardColors()[which(!standardColors() %in% omit)]
  if(length(colors) < n) n = length(colors)
  return(colors[1:n])
}

labels2colors_dt = function(vector, ...){
  return(labels2colors(vector, colorSeq = dtColors(length(vector), ...)))
}

ClusterDendrogram = function(object, assay = "RNA", group.by = "seurat_clusters"){
  Idents(object) = group.by
  DefaultAssay(object) = assay
  # object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  object = BuildClusterTree(object = object)
  return(object)
}

ComputeM1 = function(original_clusters, permuted_clusters){
  
  tableMatrix = table(original_clusters, permuted_clusters)
  
  # Proportion of type in each cluster (row-normalized)
  rnMatrix = tableMatrix/rowSums(tableMatrix)
  
  m1.scores = sapply(seq_along(unique(original_clusters)), function(index){
    proportions = rnMatrix[index,]
    score = 1-ShannonEntropy(proportions)/log2(length(proportions))
    return(score)
  })
  
  return(m1.scores)
}

ComputeM2 = function(original_clusters, permuted_clusters){
  tableMatrix = table(original_clusters, permuted_clusters)
  
  # Proportion of cluster from each type (column-normalized)
  cnMatrix = t(t(tableMatrix)/colSums(tableMatrix))
  
  m2.scores = sapply(seq_along(unique(original_clusters)), function(index){
    raw_score = sum((tableMatrix[index,] * (1-cnMatrix[index,])))/sum(tableMatrix[index,])
    score = 1-raw_score
    return(score)
  })
  
  return(m2.scores)
}

MatchClusters = function(reference, target, log.p = 10, overlap = 30, jaccard = 0){
  
  # Overlap statistics
  overlap.stats = OverlapStatistics(table(reference, target))
  overlap.signif = subset(overlap.stats, pval < 1e-10)
  
  # Throw out multi-mapping clusters
  message('Found the following multimapping clusters')
  print(subset(overlap.signif, ident1 %in% getDuplicates(overlap.signif$ident1)) %>% arrange(ident1))
  overlap.signif.one2one = subset(overlap.signif, !ident1 %in% getDuplicates(overlap.signif$ident1))
  
  # Assign new names
  transferred = convert_values(reference, key = overlap.signif.one2one[,c('ident1', 'ident2')] %>% setNames(c('old.names', 'new.names')))
  return(transferred)
}

# Given an OrthoType cluster, it finds the corresponding species clusters that match it given the thresholds 
CorrespondingCluster = function(overlap.list, cluster, log.p = 10, overlap = 30, jaccard = 0){
  if(inherits(overlap.list, 'list')){
    filtered.list = lapply(overlap.list, function(table) table[table$log.p >= log.p & 
                                                                 table$overlap >= overlap & 
                                                                 table$jaccard >= jaccard,])
    results = do.call(rbind, lapply(filtered.list, function(x) subset(x, ident1 == cluster)))
    
  } else if(inherits(overlap.list, 'data.frame')) {
    table = overlap.list
    filtered.list = table[table$log.p >= log.p & table$overlap >= overlap & table$jaccard >= jaccard,]
    results = subset(filtered.list, ident1 == cluster)
  } else {
    stop('should be a data.frame or list')
  }
  
  rownames(results) = NULL
  return(results)
}

OverlapStatistics = function(table){
  pval = HyperMatrix(table, log.p = FALSE)
  log.p = HyperMatrix(table, log.p = TRUE)
  jaccard = JSMatrix(table)
  
  # Summary table
  data=reshape2::melt(table)
  data$pval=reshape2::melt(pval)$value
  data$log.p=reshape2::melt(log.p)$value
  data$jaccard=reshape2::melt(jaccard)$value
  colnames(data)=c("ident1", "ident2", "overlap", 'pval', "log.p", "jaccard")
  return(data)
}

HyperMatrix = function(table, log.p = TRUE){
  pval = table
  for(column in 1:ncol(table)){
    for(row in 1:nrow(table)){
      pval[row,column]=phyper(table[row,column]-1,
                              sum(table[row,]),
                              sum(rowSums(table))-sum(table[row,]),
                              sum(table[,column]),
                              lower.tail=FALSE, 
                              log.p=FALSE)
    }
  }
  
  # Return negative log10 p-value
  if(log.p) pval = -log10(pval)
  
  return(pval)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

TitlePlot = function(plot, title, hjust = 0.5){
  plot + ggtitle(title) + theme(plot.title = element_text(hjust = hjust, face = 'plain'))
}

TransferMetadata = function(from, to){
  nCount_RNA = to$nCount_RNA
  nFeature_RNA = to$nFeature_RNA
  to@meta.data = from@meta.data
  to$nCount_RNA = nCount_RNA
  to$nFeature_RNA = nFeature_RNA
  
  return(to)
}

UpdateEnrichment = function(object, species, enrichment.order = c("CD90+", "CD73-", "NEUN+", "NEUN-", "CHX10+", "NEUN-\nCHX10-", "NEUN-\nCHX10+", "NEUN+\nCHX10-", "CHX10-\nCD73-\nCD133-", "NONE")){
  
  # Read in species metadata
  species_metadata = read.csv("../../Species_Objects/sample_metadata2.csv")
  
  if(species == "Human") {
    object$region = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,1])
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,2])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,3]))))
    object$donor = object$animal
  } else if(species == "Goldfish"){
    object$enrichment = "NONE"
    object$animal = object$orig.ident
  } else if(species == "Macaque") {
    # Region 
    object$region = "Periphery"
    object$region[grepl("Fovea", object$orig.ident, fixed = TRUE)] = "Fovea"
    
    # Animal
    object$animal = substr(object$tag, 0, 2)
    object$animal[startsWith(object$tag, "M1Per")] = "M5"
    object$animal[startsWith(object$tag, "M2Per")] = "M6"
    object$animal[startsWith(object$tag, "M3Per")] = "M7"
    
    # Enrichment type
    object$enrichment = "NONE"
    object$enrichment[grepl("CD90|Mixed", object$tag)] = "CD90+"
    object$enrichment[grepl("CD73-", object$tag)] = "CD73-"
    
    # object$batch = paste0(object$region, "_", object$animal, "_", object$enrichment)
    object$batch = paste0(object$region, "_", object$animal)
  } else if(species == "Marmoset") {
    object$region = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,1])
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,2])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,3]))))
  } else {
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 2)[,1])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 2)[,2]))))
  }
  
  # Make sure they are all uniform
  object$enrichment[toupper(object$enrichment) == "CD90"] = "CD90+"
  object$enrichment[toupper(object$enrichment) == "CHX10"] = "CHX10+"
  object$enrichment[toupper(object$enrichment) == "CHX10+"] = "CHX10+"
  object$enrichment[toupper(object$enrichment) == "NEUN"] = "NEUN+"
  object$enrichment[toupper(object$enrichment) == "NEUN+"] = "NEUN+"
  object$enrichment[toupper(object$enrichment) == "NEUN-"] = "NEUN-"
  object$enrichment[toupper(object$enrichment) == "CD73"] = "CD73-"
  object$enrichment[toupper(object$enrichment) == "ALL"] = "NONE"
  object$enrichment[toupper(object$enrichment) == "ISL2B"] = "ISL2B+"
  object$enrichment[toupper(object$enrichment) == "NEUN-CHX10-"] = "NEUN-\nCHX10-"
  object$enrichment[toupper(object$enrichment) == "NEUN-CHX10+"] = "NEUN-\nCHX10+"
  object$enrichment[toupper(object$enrichment) == "NEUN+CHX10-"] = "NEUN+\nCHX10-"
  
  # Order
  object$enrichment = factor(object$enrichment, levels = enrichment.order)
  
  return(object)
}

ReadAmacrineData = function(speciesList, mc.cores = 1){
  objectList = mclapply(speciesList, function(species) {
    if(species == "Mouse") {
      # return(readRDS("../../Species_Reference/MouseACref_v4.rds"))
      
      YanAC = readRDS("../../Species_Reference/YanAC_v4.rds")
      
      # Yan et al. cell depletion strategy
      # BCs and Müller glia - GFP, rods - CD73, cones - CD133
      # objectList$Mouse$enrichment = "NONE"
      YanAC$enrichment = "CHX10-\nCD73-\nCD133-"
      YanAC$seurat_clusters = YanAC$cluster_no
      
      # Assign literature types
      yan_otherdata = fread("../../data/MouseAC_other_meta.txt", sep = " ", fill = TRUE)
      YanAC$lit_type = yan_otherdata$Notes[match(YanAC$cluster_no, yan_otherdata$Cluster)]
      YanAC$lit_type[YanAC$cluster_no == 10] = "nGnG-2_CCK"
      YanAC$lit_type[YanAC$cluster_no == 63] = "PENK_SST"
      YanAC$lit_type[YanAC$lit_type == ""] = NA
      
      # Original classification
      YanAC$orig.classification = ifelse(YanAC$classification == "GABAergic", "GABA", 
                                         ifelse(YanAC$classification == "Glycinergic", "Gly", 
                                                ifelse(YanAC$classification == "Both", "Both", "nGnG")))
      
      return(YanAC)
    } else {
      readRDS(paste0("../../Species_Objects/", species, "AC_v5.rds"))
    }
  }, mc.cores = mc.cores)
  names(objectList) = speciesList
  
  return(objectList)
}

LogAvgExpr = function(object, ...){
  args = list(...)
  return(as.data.frame(log1p(AverageExpression(object, ...)[[args$assay]])))
}

GetAnnotationColors = function(vector, annotation){
  return(Colors(annotation)[match(vector, Annotation(annotation))])
}

FormatTitle = function(){
  return(theme(plot.title = element_text(hjust = 0.5, face = "plain")))
}

theme_cowplot2 = function (font_size = 14, font_family = "", line_size = 0.5, 
                           rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14) 
{
  half_line <- font_size/2
  small_size <- rel_small * font_size
  theme_grey(base_size = font_size, base_family = font_family) %+replace% 
    theme(line = element_line(color = "black", size = line_size, 
                              linetype = 1, lineend = "butt"), rect = element_rect(fill = NA, 
                                                                                   color = NA, size = line_size, linetype = 1), text = element_text(family = font_family, 
                                                                                                                                                    face = "plain", color = "black", size = font_size, 
                                                                                                                                                    hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, 
                                                                                                                                                    margin = margin(), debug = FALSE), axis.line = element_line(color = "black", 
                                                                                                                                                                                                                size = line_size, lineend = "square"), axis.line.x = NULL, 
          axis.line.y = NULL, axis.text = element_text(color = "black", 
                                                       size = small_size), axis.text.x = element_text(margin = margin(t = small_size/4), 
                                                                                                      vjust = 1), axis.text.x.top = element_text(margin = margin(b = small_size/4), 
                                                                                                                                                 vjust = 0), axis.text.y = element_text(margin = margin(r = small_size/4), 
                                                                                                                                                                                        hjust = 1), axis.text.y.right = element_text(margin = margin(l = small_size/4), 
                                                                                                                                                                                                                                     hjust = 0), axis.ticks = element_line(color = "black", 
                                                                                                                                                                                                                                                                           size = line_size), axis.ticks.length = unit(half_line/2, 
                                                                                                                                                                                                                                                                                                                       "pt"), axis.title.x = element_text(margin = margin(t = half_line/2), 
                                                                                                                                                                                                                                                                                                                                                          vjust = 1), axis.title.x.top = element_text(margin = margin(b = half_line/2), 
                                                                                                                                                                                                                                                                                                                                                                                                      vjust = 0), axis.title.y = element_text(angle = 90, 
                                                                                                                                                                                                                                                                                                                                                                                                                                              margin = margin(r = half_line/2), vjust = 1), 
          axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2), 
                                            vjust = 0), legend.background = element_blank(), 
          legend.spacing = unit(font_size, "pt"), legend.spacing.x = NULL, 
          legend.spacing.y = NULL, legend.margin = margin(0, 
                                                          0, 0, 0), legend.key = element_blank(), legend.key.size = unit(1.1 * 
                                                                                                                           font_size, "pt"), legend.key.height = NULL, legend.key.width = NULL, 
          legend.text = element_text(size = rel(rel_small)), 
          legend.text.align = NULL, legend.title = element_text(hjust = 0), 
          legend.title.align = NULL, legend.position = "right", 
          legend.direction = NULL, legend.justification = c("left", 
                                                            "center"), legend.box = NULL, legend.box.margin = margin(0, 
                                                                                                                     0, 0, 0), legend.box.background = element_blank(), 
          legend.box.spacing = unit(font_size, "pt"), panel.background = element_blank(), 
          panel.border = element_blank(), panel.grid = element_blank(), 
          panel.grid.major = NULL, panel.grid.minor = NULL, 
          panel.grid.major.x = NULL, panel.grid.major.y = NULL, 
          panel.grid.minor.x = NULL, panel.grid.minor.y = NULL, 
          panel.spacing = unit(half_line, "pt"), panel.spacing.x = NULL, 
          panel.spacing.y = NULL, panel.ontop = FALSE, strip.background = element_rect(fill = "grey80"), 
          strip.text = element_text(size = rel(rel_small), 
                                    margin = margin(half_line/2, half_line/2, half_line/2, 
                                                    half_line/2)), strip.text.x = NULL, strip.text.y = element_text(angle = -90), 
          strip.placement = "inside", strip.placement.x = NULL, 
          strip.placement.y = NULL, strip.switch.pad.grid = unit(half_line/2, 
                                                                 "pt"), strip.switch.pad.wrap = unit(half_line/2, 
                                                                                                     "pt"), plot.background = element_blank(), plot.title = element_text(
                                                                                                       size = rel(rel_large), hjust = 0.5, vjust = 1, 
                                                                                                       margin = margin(b = half_line)), plot.subtitle = element_text(size = rel(rel_small), 
                                                                                                                                                                     hjust = 0, vjust = 1, margin = margin(b = half_line)), 
          plot.caption = element_text(size = rel(rel_tiny), 
                                      hjust = 1, vjust = 1, margin = margin(t = half_line)), 
          plot.tag = element_text(face = "bold", hjust = 0, 
                                  vjust = 0.7), plot.tag.position = c(0, 1), plot.margin = margin(half_line, 
                                                                                                  half_line, half_line, half_line), complete = TRUE)
}

SmartSaveRDS = function(object, path, overwrite = FALSE){
  if(!basename(path) %in% list.files(dirname(path))){
    message("Writing file to ", path)
    saveRDS(object, path)
  } else if(overwrite){
    message("Overwriting file ", path)
    saveRDS(object, path)
  }
  else {
    warning("File already exists, skipped writing file...")
  }
}

RotatedXAxis = function(){
  return(theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
}

PrepStack = function(plt.list, nrow = NULL, ncol = NULL){
  
  # Remove x axis labels for all but last row
  plt.list[1:(length(plt.list)-ncol)] = lapply(plt.list[1:(length(plt.list)-ncol)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    return(plot)
  })
  
  # Remove y axis labels for all but first column
  if(ncol > 1){
    plt.list[1:length(plt.list) %% ncol != 1] = lapply(plt.list[1:length(plt.list) %% ncol != 1], function(plot) {
      plot = plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
      return(plot)
    })
  }
  
  return(plt.list)
}

StackedPlots = function(plt.list, height = 2, remove.titles = TRUE, no.legend = FALSE, combine.legend = TRUE, ...){
  
  # Modify plots
  # plt.list = lapply(plt.list, function(plot) {
  #   gene = plot$labels$title
  #   plot = plot + ylab(gene) + theme(axis.title.x = element_blank(), plot.title = element_blank())
  #   return(plot)
  # })
  
  # Remove axis labels for all but last plot
  plt.list[1:(length(plt.list)-1)] = lapply(plt.list[1:(length(plt.list)-1)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    return(plot)
  })
  
  # Remove legend from all except middle panel
  if(combine.legend) {
    middle.element = ceiling(length(plt.list)/2)
    plt.list[-middle.element] = lapply(plt.list[-middle.element], function(plot) plot + NoLegend()) 
  }
  
  # Remove all titles
  if(remove.titles){
    plt.list = lapply(plt.list, function(plot) {
      plot = plot + theme(plot.title = element_blank()) 
      return(plot)
    })
  }
  
  # Remove legend
  if(no.legend){
    plt.list = lapply(plt.list, function(plot) {
      plot = plot + NoLegend() 
      return(plot)
    })
  }
  
  # plt = ggarrange(plotlist = plt.list,
  #                 nrow = length(plt.list),
  #                 # align = "v",
  #                 # axis = "bt",
  #                 heights = c(rep(1, length(plt.list)-1), height),
  #                 ...)
  
  plt = patchwork::wrap_plots(plotlist = plt.list, ncol = 1)
  
  return(plt)
}

ViolinPlotFlipped = function(object, height = 2, ...){
  plt.list = VlnPlot(object, ..., combine = FALSE) 
  
  # Flip coordinates
  plt.list = lapply(plt.list, function(plot) plot + scale_y_continuous(position = "right") + coord_flip())
  
  # Modify plots
  plt.list = lapply(plt.list, function(plot) {
    gene = plot$labels$title
    plot = plot + ylab(gene) + theme(axis.title.y = element_blank(), 
                                     plot.title = element_blank(), 
                                     legend.position = "none")
    return(plot)
  })
  
  # Remove axis labels for all but last plot
  plt.list[2:(length(plt.list))] = lapply(plt.list[2:(length(plt.list))], function(plot) {
    plot = plot + theme(axis.text.y = element_blank())
    return(plot)
  })
  
  # plt = ggarrange(plotlist = plt.list, 
  #           common.legend = TRUE, 
  #           nrow = length(plt.list), 
  #           legend = 'right', 
  #           align = "v", 
  #           heights = c(rep(1, length(plt.list)-1), height))
  
  plt = patchwork::wrap_plots(plotlist = plt.list, nrow = 1)
  return(plt)
}

ViolinPlot = function(object, height = 2, ...){
  plt.list = VlnPlot(object, ..., combine = FALSE)
  
  # Modify plots
  plt.list = lapply(plt.list, function(plot) {
    gene = plot$labels$title
    plot = plot + ylab(gene) + theme(axis.title.x = element_blank(), plot.title = element_blank())
    return(plot)
  })
  
  # Remove axis labels for all but last plot
  plt.list[1:(length(plt.list)-1)] = lapply(plt.list[1:(length(plt.list)-1)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank())
    return(plot)
  })
  
  plt = ggarrange(plotlist = plt.list, 
                  common.legend = TRUE, 
                  nrow = length(plt.list), 
                  legend = 'right', 
                  align = "v", 
                  heights = c(rep(1, length(plt.list)-1), height))
  return(plt)
}

GetSpecies = function(object, field = "species"){
  return(unique(object@meta.data[[field]]))
}

# Annotate clusters using a dataframe with old and new names
AnnotateClusters = function(object, key, from = "seurat_clusters", to = "annotated", factor = FALSE){
  new.names = key$new.names[match(as.character(object@meta.data[,from]), as.character(key$old.names))]
  if(factor) return(factor(new.names, levels = sort(unique(key$new.names)))) else return(new.names)
}

theme_umap = function(plt, remove.axes = FALSE, title = NULL, rename.axes = TRUE){
  
  # Correct xlab ylab
  if(rename.axes) plt = plt + labs(x = 'UMAP 1', y = 'UMAP 2')
  
  if(remove.axes) {
    plt2 = plt + theme_void()
  } else{
    plt2 = plt + theme_dario() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank())
  }
  
  # Center title 
  plt2 = plt2 + 
    {if(!is.null(title)) ggtitle(title)}+ 
    theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.title=element_blank())
  
  return(plt2)
}

TransferSACLabels = function(object){
  AC_model <- TrainXGBoost(PengAC, object, train.clusters = "lit_type")
  # saveRDS(AC_model, "../../data/train_Peng_test_Rat.rds")
  # train_Peng_test_Rat <- BuildConfusionMatrix(object, PengAC, model = AC_model)
  object$xgb.labels = PredictLabels(object, model = AC_model, scale.by.model = TRUE)
  object$seurat_clusters = factor(object$seurat_clusters)
  print(JSHeatmap(JSMatrix(table(object$xgb.labels, object$seurat_clusters)), heatmap = FALSE))
  # print(plotConfusionMatrix(table(as.character(object$xgb.labels), object$seurat_clusters)))
  return(object)
}

plotSACMarkers = function(object, group.by = "seurat_clusters", feature.plot = FALSE){
  if(feature.plot){
    plot_grid(
      {if(any(ON_SAC_markers %in% rownames(object))) FeaturePlot(object, features = ON_SAC_markers, ncol = length(ON_SAC_markers))}, 
      {if(any(OFF_SAC_markers %in% rownames(object))) FeaturePlot(object, features = OFF_SAC_markers, ncol = length(OFF_SAC_markers))}, 
      ncol = 1, nrow = 2)
  } else {
    plot_grid(
      {if(any(ON_SAC_markers %in% rownames(object))) VlnPlot(object, group.by = group.by, features = ON_SAC_markers, ncol = length(ON_SAC_markers))}, 
      plot_grid(
        {if(any(OFF_SAC_markers %in% rownames(object))) VlnPlot(object, group.by = group.by, features = OFF_SAC_markers, ncol = length(OFF_SAC_markers))}, 
        VlnPlot(object, features = c("ON.SAC.score", "OFF.SAC.score")),
        ncol = 2, rel_widths = c(3,2)
      ),
      ncol = 1, nrow = 2)
  }
}

plotXGBLabels = function(object){
  plot_grid(
    plot_grid(
      DimPlot(object, group.by = "seurat_clusters", label = TRUE) + NoLegend(), 
      DimPlot(object, group.by = "xgb.labels", label = TRUE) + NoLegend()
    ),
    plotSACMarkers(object, group.by = "xgb.labels"), 
    ncol = 1, nrow = 2)
}

SubclusterSACs = function(object, group.by = "animal", resolution = 0.2, harmonize = TRUE){
  if(harmonize == TRUE){
    object = Harmonize(object, batch = group.by, cluster_resolution = resolution, show.plots = FALSE)
  } else {
    object = ClusterSeurat(object, cluster_resolution = resolution)
  }
  print(BrowseSeurat(object, batch = group.by))
  object = TransferSACLabels(object)
  print(plotSACMarkers(object))
  print(plotSACMarkers(object, feature.plot = TRUE))
  return(object)
}

BarPlot = function(data, x, y, color = NULL){
  p = ggplot(data, aes(x=eval(parse(text = x)), y=eval(parse(text = y)))) + 
    geom_col(color = "black") +
    xlab(x)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    ylab(y)+
    theme_cowplot()+
    RotatedAxis()
  
  return(p)
}

stackedBarGraph2 = function(data, x, y, proportion = TRUE, position = position_stack(), label = FALSE, label.threshold = 0.05){
  
  if(inherits(data, "Seurat")) data = data@meta.data
  
  # data[,x] = as.character(data[,x])
  # data[,y] = as.character(data[,y])
  cross.tabulation = table(data[,x], data[,y])
  
  if(proportion) cross.tabulation = cross.tabulation/rowSums(cross.tabulation)
  # ggbarplot(reshape2::melt(cross.tabulation), x = "Var1", y = "value", fill = "Var2", position = position, label = reshape2::melt(cross.tabulation)[["Var2"]], lab.col = "black", lab.pos = "in") + RotatedAxis() +
  #   xlab(x)+ ylab(y)
  melted = reshape2::melt(cross.tabulation)
  melted$Var2 = as.character(melted$Var2)
  melted$label = melted$Var2
  melted$label[melted$value < label.threshold] = NA
  
  p = ggplot(melted, aes(fill=Var2, y=value, x=Var1, label = label)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    theme_minimal()+
    # xlab(x)+
    labs(fill=y)+
    {if(label) geom_text(size = 3, position = position_stack(vjust = 0.5))}+
    # scale_x_continuous(breaks = seq(0, nClusters(seurat)-1, 1))+
    # geom_text_repel()+
    # scale_fill_manual(feature.1)+
    # scale_x_discrete()
    scale_y_continuous(expand = c(0,0))+
    ylab('Proportion')+
    theme_cowplot()+
    # theme_dario()+
    RotatedAxis() + 
    theme(axis.title.x = element_blank())
  # theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))
  return(p)
}

CorrelationHeatmap = function(object, 
                              assay = "RNA", 
                              group.by = "seurat_clusters", 
                              annotate.by = NULL, 
                              cluster_rows = TRUE, 
                              features = NULL, 
                              title = NULL, 
                              label = FALSE, 
                              cols = c("white", "red"), 
                              annotation_cols = NULL, 
                              return.dendrogram = FALSE, 
                              ...){
  
  if(!is.null(features)){
    avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features)[[assay]]))
  } else {
    message("Using top 2000 variable features!")
    object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
    avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = object@assays[[assay]]@var.features)[[assay]]))
  }
  print(dim(avg_exp))
  type_cor = cor(avg_exp)
  
  if(!is.null(annotate.by)){
    metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by)
    annotation = metadata[match(colnames(type_cor), metadata[[group.by]]), annotate.by]
    # annotation_cols2 = c(annotation_cols[names(annotation_cols) %in% annotation])
    if(is.null(annotation_cols)){
      ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left")
    } else{
      ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left", col = list(anno = annotation_cols))
    }
  } else {
    ha = NULL
  }
  
  if(inherits(cluster_rows, 'dendrogram')){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, 
                 cluster_rows  = as.hclust(cluster_rows),
                 cluster_columns = as.hclust(cluster_rows),
                 border_gp = gpar(col = "black", lty = 1),
                 #show_row_names = TRUE, show_column_names = TRUE, row_dend_side = "left", 
                 show_column_dend = TRUE, 
                 column_title = title,
                 cell_fun = if(label) function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", type_cor[i, j]), x, y, gp = gpar(fontsize = 10))} else NULL,
                 heatmap_legend_param = list(title = "Pearson R"),
                 top_annotation = ha, 
                 ...)
    ht = draw(ht)
  } else if(cluster_rows == TRUE){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, #colorRamp2(seq(0,1,0.1), rev(rainbow(11))), #colorRamp2(c(-1,0,1), c("blue", "white", "red")),
                 cluster_rows = hclust(as.dist(1 - type_cor)),
                 cluster_columns = hclust(as.dist(1 - type_cor)), 
                 border_gp = gpar(col = "black", lty = 1),
                 #show_row_names = TRUE, show_column_names = TRUE, row_dend_side = "left", 
                 show_column_dend = TRUE, 
                 column_title = title,
                 cell_fun = if(label) function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", type_cor[i, j]), x, y, gp = gpar(fontsize = 10))} else NULL,
                 heatmap_legend_param = list(title = "Pearson R"),
                 top_annotation = ha, 
                 ...)
    ht = draw(ht)
  } else if(cluster_rows == FALSE){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, #colorRamp2(seq(0,1,0.1), rev(rainbow(11))),
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE, 
                 show_row_names = TRUE, 
                 show_column_names = TRUE, 
                 row_dend_side = "left", 
                 show_column_dend = TRUE, 
                 column_title = title,
                 heatmap_legend_param = list(title = "Pearson R"),
                 top_annotation = ha, ...)
    ht = draw(ht)
  } else {
    stop('could not parse cluster_rows')
  }
  
  if(return.dendrogram) row_dend(ht)
}

EvaluateModel = function(table){
  TP = table[1,1]
  TN = table[2,2]
  FP = table[1,2]
  FN = table[2,1]
  accuracy = (TP+TN)/sum(rowSums(table))
  sensitivity = (TP)/(TP + FN) # true positive rate
  specificity = (TN)/(TN + FP) # true negative rate
  precision = TP/(TP+FP) # proportion of predicted positives that are truly positive
  recall = TP/(TP+FN) # same as sensitivity
  accuracy = (TP+TN)/sum(rowSums(table))
  message(paste0(capture.output(table), collapse = '\n'))
  message("accuracy: ", accuracy, "\ntrue positive rate (AKA sensitivity or recall): ", sensitivity, "\nfalse positive rate: ", 1-specificity, '\nprecision: ', precision)
}

ScatterHistogram = function(object, x, y, color = NULL, offset = -0.2){
  
  if(inherits(object, "Seurat")){
    data = object@meta.data
  } else {
    data = object
  }
  
  p1 = gghistogram(data, x = x, xlab = FALSE) + theme(axis.text.x = element_blank())
  p2 = NULL
  p3 = ggscatter(data, x, y, color = color) + theme(legend.position = "bottom")
  p4 = gghistogram(data, x = y, ylab = FALSE) + coord_flip() + RotatedAxis() + theme(axis.text.y = element_blank())
  
  plot_grid(
    p1, NULL, NULL, 
    NULL, NULL, NULL, 
    p3, NULL, p4,
    ncol = 3, 
    nrow = 3, 
    rel_widths = c(2, offset, 1), 
    rel_heights = c(1, offset, 2), 
    align = "hv", 
    axis = "btlr")
  
  # patchworked
  # (p1 | plot_spacer()) / (p3 | p4)
}

MixtureHistogram = function(df, variable, model, intersects = NULL, color.1 = "magenta", color.2 = "blue", legend.name = NULL, xlab = NULL, ylab = "Density", plot.points = TRUE, logY = FALSE){
  ggplot(df, aes(x=eval(parse(text = variable)))) +
    geom_histogram(bins = 20, aes(y=..density..), colour="black", fill = "white", boundary=0) + 
    # geom_density(aes(y=..density..)) +
    geom_function(fun = function(x) model$lambda[1]*dnorm(x, mean = model$mu[1], sd = model$sigma[1]), linewidth = 1, linetype = "solid", color = color.1)+
    geom_function(fun = function(x) model$lambda[2]*dnorm(x, mean = model$mu[2], sd = model$sigma[2]), linewidth = 1, linetype = "solid", color = color.2)+
    geom_vline(xintercept = c(intersects), linetype = "dashed", color = "grey")+
    scale_x_continuous(breaks = seq(0,1, by = 0.25), expand = c(0,0))+
    {if(plot.points) geom_point(aes(fill = eval(parse(text = legend.name)), y = 0), shape = 21, color = "black", size = 3)}+
    scale_fill_gradient(name = legend.name, limits = c(0,1), low = color.1, high = color.2)+
    {if(logY) scale_y_continuous(trans = "log10")}+
    ylab(ylab)+
    xlab(xlab)+
    theme_dario()
}

ComputePosteriors = function(xi, model){
  sigma1 = model$sigma[1]
  sigma2 = model$sigma[2]
  mu1 = model$mu[1]
  mu2 = model$mu[2]
  lambda1 = model$lambda[1]
  lambda2 = model$lambda[2]
  
  p_xi_1 = dnorm(xi, mean = mu1, sd = sigma1)
  p_1 = lambda1
  
  p_xi_2 = dnorm(xi, mean = mu2, sd = sigma2)
  p_2 = lambda2
  
  # Bayes Theorem
  results.df = data.frame(p_1_xi = (p_xi_1 * p_1)/(p_xi_2 * p_2 + p_xi_1 * p_1), 
                          p_2_xi = (p_xi_2 * p_2)/(p_xi_2 * p_2 + p_xi_1 * p_1))
  results.df$check = results.df$p_1_xi + results.df$p_2_xi
  
  return(results.df)
}

FindCutpoint = function(model){
  sigma1 = model$sigma[1]
  sigma2 = model$sigma[2]
  mu1 = model$mu[1]
  mu2 = model$mu[2]
  alpha = model$lambda[1]
  beta = model$lambda[2]
  
  # Intersection point of two Gaussians is where probability flips (adapted from https://stats.stackexchange.com/questions/311592/how-to-find-the-point-where-two-normal-distributions-intersect)
  A = (-1/sigma1^2+1/sigma2^2)
  B = 2*(-mu2/sigma2^2 + mu1/sigma1^2)
  C = (mu2^2/sigma2^2) - (mu1^2/sigma1^2) + log(sigma2^2/sigma1^2) + 2*log(alpha/beta)
  discriminant = B^2 - 4*A*C
  intersect1 = (-B + sqrt(discriminant))/(2*A)
  intersect2 = (-B - sqrt(discriminant))/(2*A)
  intersect = c(intersect1, intersect2)[which(c(intersect1, intersect2) > min(c(mu1,mu2)) & c(intersect1, intersect2) < max(c(mu1,mu2)))]
  return(intersect)
}

UpdateLitTypes = function(object, field = "lit_type"){
  previous_names = object@meta.data[[field]]
  object@meta.data[[field]][previous_names == "AII"] = "A2"
  object@meta.data[[field]][previous_names == "CA-I"] = "CA1"
  object@meta.data[[field]][previous_names == "CA-II"] = "CA2"
  return(object)
}

ConservedEnrichedHeatmap2 = function(object1, object2, object3, geneList, color, ...){
  
  p1 = DoHeatmap(object1, features = unlist(geneList), ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  p2 = DoHeatmap(object2, features = unlist(geneList), ...) + 
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  p3 = DoHeatmap(object3, features = unlist(geneList), ...) + 
    NoLegend() +
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  panel = plot_grid(p1, NULL, p2, NULL, p3, 
                    # rel_widths = c(3, -0.05, 3, -0.05, 3),
                    rel_widths = c(3, 0, 3, 0, 3),
                    ncol = 5) 
  
  return(panel)
}

ConservedEnrichedHeatmap = function(object, sharedList, primateList, rodentList, laurasiaList, ...){
  # Shared
  s1 = DoHeatmap(object, features = (unlist(sharedList)),  ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "magenta"))
  # coord_flip()
  
  # Species-specific 
  p1 = DoHeatmap(object, features = unlist(primateList),  ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "blue"))
  # coord_flip()
  
  p2 = DoHeatmap(object, features = unlist(rodentList),  ...) + 
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradientn(colors = c("white", "white", "turquoise"))
  # coord_flip()
  
  p3 = DoHeatmap(object, features = unlist(laurasiaList),  ...) + 
    NoLegend() +
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "orange"))
  # coord_flip()
  
  panel = plot_grid(s1,p1,p2,p3, 
                    ncol = 4) 
  
  return(panel)
}

theme_dario = function(margin = 10){
  dario.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_rect(fill = NA, colour = "black", linewidth=1),
                      panel.background = element_blank(), 
                      plot.margin = margin(t = margin, r = margin, b = margin, l = margin),
                      axis.line = element_blank(), 
                      # panel.background = element_rect(fill = 'white', colour = "black", linewidth=1), 
                      axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), 
                      plot.title = element_text(hjust = 0.5)) 
  return(dario.theme)
}

PlotCelltypeProportions = function(objectList, proportions, celltype, color = "orange", reference.line){
  
  SAC.clusters = lapply(objectList, function(object) {
    metadata = Metadata(object, feature.1 = "seurat_clusters", feature.2 = "lit_type")
    return(metadata$seurat_clusters[which(metadata$lit_type == celltype)])
  })
  
  SAC.clusters = SAC.clusters[match(names(proportions), names(SAC.clusters))]
  
  SAC.props = data.frame(species = factor(c(names(SAC.clusters)), levels = phylogenetic_order), 
                         proportion = 100*c(unlist(sapply(seq_along(proportions), function(i) sum(proportions[[i]][SAC.clusters[[i]] ])))), 
                         group = c(rep("Test", length(objectList))))
  
  ggplot(SAC.props, aes(x = species, y = proportion, fill = group)) + 
    ylab(paste0(celltype, " proportion (%)"))+
    geom_col(color = "black")+
    geom_hline(yintercept = reference.line, linetype = "dashed")+
    scale_fill_manual(values = c(color))+
    scale_y_continuous(expand = c(0,0))+
    theme_cowplot()+
    theme(legend.position = "none")+
    RotatedAxis() 
  # coord_flip()
}

# Downsample seurat by variable identity class
DownsampleSeurat = function(object, group.by, size, seed = 12345){
  # previous_ident = Idents(object)
  Idents(object) = group.by
  object = subset(object, downsample = size, seed = seed)
  # object = subset(object, cells = WhichCells(seurat, downsample = size, seed = seed))
  # Idents(object) = previous_ident
  return(object)
}

ComputeIntegrationLISI  = function(object, group.by = "species", plot = TRUE, subsample.types = TRUE, nPCs = 20, method = 'seurat'){
  
  if(subsample.types) {
    # Downsample to minimum category size
    downsample = min(table(object[[group.by]]))
    object = DownsampleSeurat(object, group.by, downsample)
  }
  
  if(method == 'seurat'){
    iLISI_observed <- compute_lisi(object@reductions$pca@cell.embeddings[,1:nPCs], 
                                   object@meta.data, c(group.by))
  } else if(method == 'harmony') {
    iLISI_observed <- compute_lisi(object@reductions$harmony@cell.embeddings[,1:nPCs], 
                                   object@meta.data, c(group.by))
  } else {
    stop('method should be seurat or harmony')
  }
  
  iLISI_expected = SimpsonIndex(table(object@meta.data[[group.by]]), invert = TRUE)
  
  # Normalized iLISI
  object$normalized.iLISI = (iLISI_observed[[group.by]]-1)/(iLISI_expected-1)
  
  if(plot){
    print(plot_grid(ClusterBatchPlot(object, batch = group.by, shuffle = TRUE), 
                    FeaturePlot(object, features = "normalized.iLISI"), 
                    VlnPlot(object, features = "normalized.iLISI", pt.size = 0)+NoLegend(),
                    VlnPlot(object, features = "normalized.iLISI", group.by = group.by, pt.size = 0)+NoLegend(),
                    ncol = 2, nrow = 2))
  }
  
  return(object)
}

ComputeLittypeLISI = function(object, types.to.consider = c("A2", "VG3", "SAC", "A17"), subset.types = FALSE, subsample.types = TRUE, nPCs = 20){
  knownAC = object
  knownAC$species_littype[!knownAC$species_littype %in% types.to.consider] = "Other"
  
  if(subset.types) knownAC = subset(knownAC, species_littype %in% types.to.consider)
  if(subsample.types) {
    # Downsample to minimum category size
    downsample = min(table(knownAC[["species_littype"]]))
    knownAC = DownsampleSeurat(knownAC, "species_littype", downsample)
  }
  
  cLISI_observed <- compute_lisi(knownAC@reductions$pca@cell.embeddings[,1:nPCs], 
                                 knownAC@meta.data, c('species_littype'))
  
  cLISI_expected = SimpsonIndex(table(knownAC$species_littype), invert = TRUE) # length(unique(knownAC$species_littype)) 
  
  # Normalized cLISI
  knownAC$normalized.cLISI = (cLISI_observed$species_littype-1)/(cLISI_expected-1)
  print(plot_grid(DimPlot(knownAC, group.by = "species_littype", label = TRUE, shuffle = TRUE) + NoLegend(), 
                  ClusterBatchPlot(knownAC, batch = "species", group.by = "species_littype", shuffle = TRUE), 
                  ClusterFeaturePlot(knownAC, features = "normalized.cLISI", group.by = "species_littype"), 
                  plot_grid(VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species_littype", pt.size = 0, log = F) + NoLegend(),
                            VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species", pt.size = 0, log = F) + NoLegend(), 
                            ncol = 2, align = "h", rel_widths = c(1,2)), 
                  ncol = 2, nrow = 2, rel_widths = c(1,1.2)))
  return(knownAC)
}

# ComputeLittypeLISITest = function(object, types.to.consider = c("A2", "VG3", "SAC", "A17")){
#   knownAC = object
#   knownAC$species_littype = sample(types.to.consider, ncol(knownAC), replace = TRUE)
#   
#   cLISI_observed <- compute_lisi(knownAC@reductions$umap@cell.embeddings, 
#                                  knownAC@meta.data, c('species_littype'))
#   
#   cLISI_expected = SimpsonIndex(table(knownAC$species_littype), invert = TRUE) # length(unique(knownAC$species_littype)) 
#   
#   # Normalized cLISI (no normalization for now)
#   knownAC$normalized.cLISI = cLISI_observed$species_littype # (cLISI_observed$species_littype-1)/(species_littype-1)
#   print(plot_grid(DimPlot(knownAC, group.by = "species_littype", label = TRUE, shuffle = TRUE) + NoLegend(), 
#                   ClusterBatchPlot(knownAC, batch = "species", group.by = "species_littype", shuffle = TRUE), 
#                   ClusterFeaturePlot(knownAC, features = "normalized.cLISI", group.by = "species_littype"), 
#                   plot_grid(VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species_littype", pt.size = 0, log = F) + NoLegend(),
#                             VlnPlot(knownAC, features = "normalized.cLISI", group.by = "species", pt.size = 0, log = F) + NoLegend(), ncol = 2, align = "h"), 
#                   ncol = 2, nrow = 2))
#   return(knownAC)
# }

SubsetSeuratGenes = function(object, features){
  counts = object@assays$RNA@counts
  newObject = CreateSeuratObject(counts[rownames(counts) %in% features,])
  nCount_RNA = newObject$nCount_RNA
  nFeature_RNA = newObject$nFeature_RNA
  newObject@meta.data = object@meta.data
  newObject$nCount_RNA = nCount_RNA
  newObject$nFeature_RNA = nFeature_RNA
  # newObject@meta.data[, !colnames(newObject@meta.data) %in% c("nCount_RNA", "nFeature_RNA")] = object@meta.data[, !colnames(newObject@meta.data) %in% c("nCount_RNA", "nFeature_RNA")]
  return(newObject)
}

MultigeneDotPlot = function(object, genes, ...){
  StackedPlots(plt.list = lapply(genes, function(gene) geneDotPlotFast(object, gene, ...)), remove.titles = FALSE, no.legend = FALSE)
}

geneDotPlotFast <- function(object, gene, group.by = "type", scale.within.species = TRUE, mini = FALSE, col.low = 'lightgrey', col.high = '#584B9FFF', mc.cores = 20){
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = mclapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by)$data, mc.cores = mc.cores)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = paste0(object$species, "_", object@meta.data[[group.by]])
    dat=DotPlot(object, features=gene)$data
    # metadata = Metadata(object, feature.1 = "species", feature.2 = group.by)
    dat$species = str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT = factor(str_split_fixed(dat$id, "_", 2)[,2], levels = levels(object@meta.data[[group.by]]))
  }
  
  # Plot bigger values last
  dat = dat %>% arrange(pct.exp)
  
  # Plot
  plt = ggplot(dat, aes(y=factor(species, levels = rev(levels(object$species))), x = OT, color=avg.exp.scaled, size=pct.exp))+
    geom_point()+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle(gene)+
    theme_classic()+ 
    coord_cartesian(clip = 'off')+
    guides(color = guide_colorbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed"))+
    theme(plot.title=element_text(hjust=0.5, face = 'italic'), 
          axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black", angle=45, hjust=1))+
    scale_color_gradient(low=col.low, high=col.high)
  
  if(mini){
    plt = plt + 
      ggtitle(gene)+
      theme(axis.text.x = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks = element_blank(), 
            axis.title = element_blank(),
            plot.title=element_text(hjust=0.5, face = 'italic'),
            legend.position = 'none')
  }
  
  return(plt)
}

geneDotPlot <- function(object, gene, group.by = "type", scale.within.species = TRUE){
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = lapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by)$data)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = paste0(object$species, "_", object@meta.data[[group.by]])
    dat=DotPlot(object, features=gene)$data
    # metadata = Metadata(object, feature.1 = "species", feature.2 = group.by)
    dat$species = str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT = factor(str_split_fixed(dat$id, "_", 2)[,2], levels = levels(object@meta.data[[group.by]]))
  }
  
  # Plot
  return(ggplot(dat, aes(y=factor(species, levels = rev(unique(species))), x = OT, color=avg.exp.scaled, size=pct.exp))+
           geom_point()+
           xlab("OT")+
           ylab("Species")+
           ggtitle(gene)+
           theme_classic()+ 
           theme(plot.title=element_text(hjust=0.5), 
                 axis.text.x=element_text(colour="black", angle=45, hjust=1))+
           scale_color_gradient(low="grey", high="blue"))
}

# We don't want scaling across all celltypes and across all species, we want scaling within species across celltypes
# Old geneSpeciesPlot:                                                      Vsx2180  13.5364298    87.676768           Vsx2                    Pig_RBC     0.367316626 
# DotPlot(object, ident = "Pig", features="Vsx2", group.by = "OT")$data: Vsx2     13.536430     87.67677            Vsx2                    RBC         2.50000000
OTDotPlot <- function(object, gene, celltype, scale.within.species = TRUE){
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = lapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by)$data)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = object$species_OT
    dat=DotPlot(object, features=gene)$data
    dat$species=stringr::str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT=stringr::str_split_fixed(dat$id, "_", 2)[,2]
  }
  
  dat=dat[dat$OT == as.character(celltype),]
  
  # Plot
  return(ggplot(dat, aes(y=factor(species, levels = rev(seurat_datasets)), x=features.plot, color=avg.exp.scaled, size=pct.exp))+
           geom_point()+
           xlab("Gene")+
           ylab("Species")+
           ggtitle(celltype)+
           theme_classic()+ 
           theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(colour="black", angle=45, hjust=1))+
           scale_color_gradient(low="grey", high="blue"))
  #return(dat)
}

ComputeCelltypeLISI = function(object, batch = "animal", nPCs = 20, harmony = FALSE){
  obj.list = SplitObject(object, split.by = batch)
  names(obj.list) = NULL
  cLISI.list = lapply(obj.list, function(each) {
    
    # Using first 20 PCs
    pcs = if(harmony) object@reductions$harmony@cell.embeddings[,1:nPCs] else object@reductions$pca@cell.embeddings[,1:nPCs]
    
    # Subset embedding 
    subset_pcs = pcs[match(Cells(each), rownames(pcs)),]
    
    # Compute cLISI
    batch.cLISI = compute_lisi(subset_pcs, each@meta.data, c("batch_clusters"))
    return(batch.cLISI)
  })
  cLISI.df = do.call(rbind, cLISI.list)
  return(cLISI.df[match(Cells(object), rownames(cLISI.df)), , drop = FALSE])
}

ClusterBatches = function(object, batch = "animal", resolution = 1.5){
  obj.list = SplitObject(object, split.by = batch)
  obj.list = lapply(obj.list, function(each) ClusterSeurat(each, cluster_resolution = resolution, do.umap = FALSE))
  names(obj.list) = NULL
  metadata = do.call(rbind, lapply(obj.list, function(each) each@meta.data))
  metadata$batch_clusters = paste0(metadata[[batch]], "_", metadata[["seurat_clusters"]])
  return(metadata$batch_clusters[match(Cells(object), rownames(metadata))])
}

getProportions = function(object, enrichment.group = c("ALL", "NONE", "CD73", "CD73-")){
  if(length(intersect(object$enrichment, c("ALL", "NONE", "CD73", "CD73-"))) == 0){
    message("No cells found!")
    return(list(NULL, 0))
  }
  object = subset(object, enrichment %in% enrichment.group)
  message("Subsetting to un-enriched samples: ", length(Cells(object)), " cells left")
  props = table(object$seurat_clusters)/sum(table(object$seurat_clusters))
  return(list(props, length(Cells(object))))
}

# Assigns each cluster in object 
AssignAnnotations = function(object, annotation, use = 'seurat_clusters', as = "lit_type"){
  object@meta.data[[as]] = NA
  for(index in seq_along(annotation)){
    object@meta.data[[as]][object@meta.data[[use]] %in% annotation[[index]] ] = names(annotation)[[index]]
  }
  
  # if(factor) object@meta.data[[as]] = factor(object@meta.data[[as]])
  
  return(object)
}

PlotOrthoFeature = function(object, gene, order = FALSE){
  new_genename = paste0(gene, ".")
  object[[new_genename]] = BigMatrix[gene,]
  FeaturePlot(object, feature = new_genename, raster = FALSE, order = order)
}

ExtractBarcodes = function(object, return.barcodes = FALSE, return.aliases = FALSE){
  if(all(endsWith(Cells(object), "-1"))){
    object_barcodes = gsub("-1", "", substrRight(Cells(object), 18))
    
    if(return.barcodes) return(object_barcodes)
    if(return.aliases) {
      # aliases = sapply(seq_along(object_barcodes), function(i) str_split_fixed(Cells(object)[i], object_barcodes[i], 2)[,1])
      aliases = substr(Cells(object), 1, nchar(Cells(object))-19)
      return(aliases)
    }
    
    # Remove ambiguous barcodes
    object_duplicates = names(table(object_barcodes))[table(object_barcodes) > 1]
    object = object[,!object_barcodes %in% object_duplicates]
    object = RenameCells(object, new.names = gsub("-1", "", substrRight(Cells(object), 18)))
    
  } else if(all(endsWith(Cells(object), "x"))){
    object_barcodes = gsub("x", "", str_split_fixed(Cells(object), ":", 2)[,2])
    
    if(return.barcodes) return(object_barcodes)
    if(return.aliases) {
      # aliases = sapply(seq_along(object_barcodes), function(i) str_split_fixed(Cells(object)[i], object_barcodes[i], 2)[,1])
      aliases = str_split_fixed(Cells(object), ":", 2)[,1]
      return(aliases)
    }
    
    # Remove ambiguous barcodes
    object_duplicates = names(table(object_barcodes))[table(object_barcodes) > 1]
    object = object[,!object_barcodes %in% object_duplicates]
    object = RenameCells(object, new.names = gsub("x", "", str_split_fixed(Cells(object), ":", 2)[,2]))
    
  } else {
    stop("Cell names did not match criteria!")
  }
  
  return(object)
}

MapBack = function(object1, object2, group.by = "seurat_clusters", convert.barcodes = FALSE, match.barcodes = FALSE, mc.cores = 1){
  
  if(convert.barcodes){
    object1 = ExtractBarcodes(object1)
    object2 = ExtractBarcodes(object2)
  }
  
  if(match.barcodes){
    # Search over all barcode types in object1 and see if they match all barcode types in object2
    object1$alias = ExtractBarcodes(object1, return.aliases = TRUE)
    object1$barcode = ExtractBarcodes(object1, return.barcodes  = TRUE)
    object2$alias = ExtractBarcodes(object2, return.aliases = TRUE)
    object2$barcode = ExtractBarcodes(object2, return.barcodes  = TRUE)

    # Compare barcodes
    intersects = do.call(rbind, mclapply(unique(object1$alias), function(i) {
      lapply(unique(object2$alias), function(j) {
        length(intersect(subset(object1, alias == i)$barcode, subset(object2, alias == j)$barcode))
      })
    }, mc.cores = mc.cores))
    
    rownames(intersects) = unique(object1$alias)
    colnames(intersects) = unique(object2$alias)
    
    # Find correspondence
    if(length(unique(object1$alias)) > length(unique(object2$alias))){
      # column max
      key = data.frame(object1 = rownames(intersects)[apply(intersects, 2, which.max)], 
                       object2 = colnames(intersects))
    } else {
      # row max
      key = data.frame(object1 = rownames(intersects), 
                       object2 = colnames(intersects)[apply(intersects, 1, which.max)])
    }
    
    print(key)
    
    # Match the barcodes
    object2 = RenameCells(object2, new.names = paste0(convert_values(object2$alias, key %>% setNames(c('new.names', 'old.names'))), '_', object2$barcode, '-1'))
    
    # Ensure same format in first object
    object1 = RenameCells(object1, new.names = paste0(object1$alias, '_', object1$barcode, '-1'))
  }
  
  object1$mapped.back = object2@meta.data[match(colnames(object1), colnames(object2)), group.by]
  return(object1)
}

SmartMerge = function(matrixList, by){
  
  # Find common
  common_genes = Reduce(intersect, lapply(matrixList, function(x) x[[by]]))
  
  # Subset each item to common
  matrixList = lapply(matrixList, function(item) item[item[[by]] %in% common_genes,])
  
  # Initialize merged matrix
  first = matrixList[[1]]
  merge = first[match(common_genes, first[[by]]),]
  for(i in 2:length(matrixList)){
    current = matrixList[[i]]
    new = current[match(common_genes, current[[by]]),]
    merge = cbind(merge, new) 
  }
  
  return(merge)
}

BigOrthologMatrix2 = function(species, orthologous_genes, ortho){
  
  message("Working on ", species, "!")
  
  # Read in species data
  if(species == "Mouse") {
    # speciesAC = readRDS("../../Species_Reference/MouseACref_v4.rds")
    # speciesAC = readRDS("../../../storage/HahnObjects/BCs/MouseBC_int_ann_v3.rds")
    speciesAC = readRDS("../../Species_Objects/Mouse_initial.rds")
  } else {
    speciesAC = readRDS(paste0("../../Species_Objects/", species, "_initial.rds"))
  }
  
  # Remove cells that are not in orthotype object
  cells.use = Cells(ortho)[Cells(ortho) %in% Cells(speciesAC)]
  message("Found ", length(cells.use), " cells from ortho object")
  speciesAC = speciesAC[,cells.use]
  
  # Revert gene names to original species names
  if(!is.null(speciesAC@misc$orig.features)) {
    rownames(speciesAC@assays$RNA@data) = speciesAC@misc$orig.features
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(!species %in% c("Human", "TreeShrew", "Rhabdomys")) {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, "../../Orthology/martMergeRefHuman.txt", species))
  } else if(species == "Rhabdomys") {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, "../../Orthology/martMergeRefHuman.txt", "Mouse"))
  }
  
  # Get count matrix
  count_matrix = speciesAC@assays$RNA@data
  
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthologous_genes),]
  message("Found ", nrow(subset_count_matrix), " genes out of ", length(orthologous_genes), " total orthologs. Efficiency = ", nrow(subset_count_matrix)/length(orthologous_genes))
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  ortholog_count_matrix = subset_count_matrix[match(orthologous_genes, rownames(subset_count_matrix)),]
  rownames(ortholog_count_matrix) = orthologous_genes
  # ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  # ortholog_count_matrix[1:10,1:5]
  
  return(ortholog_count_matrix)
}

BigOrthologMatrix = function(species, orthologous_genes, mart_file = '../../Orthology/martMergeRefHuman_modified.rds'){
  
  message("Working on ", species, "!")
  # Read in species data
  if(species == "Mouse") {
    speciesAC = readRDS("../../Species_Reference/MouseACref_v4.rds")
    # speciesAC = readRDS("../../Species_Reference/YanAC_v4.rds")
  } else {
    speciesAC = readRDS(paste0("../../Species_Objects/", species, "AC_v6.rds"))
  }
  
  # Revert names to original species names
  if(!is.null(speciesAC@misc$orig.features)) {
    rownames(speciesAC@assays$RNA@counts) = speciesAC@misc$orig.features # must be the counts since this is what convert ConvertGeneSymbols uses
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(!species %in% c("Human", "TreeShrew", "Rhabdomys")) {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, mart_file, species))
  } else if(species == "Rhabdomys") {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, mart_file, "Mouse"))
  }
  
  # Get count matrix
  count_matrix = speciesAC@assays$RNA@data
  
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthologous_genes),]
  message("Found ", nrow(subset_count_matrix), " genes out of ", length(orthologous_genes), " total orthologs. Efficiency = ", nrow(subset_count_matrix)/length(orthologous_genes))
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  ortholog_count_matrix = subset_count_matrix[match(orthologous_genes, rownames(subset_count_matrix)),]
  rownames(ortholog_count_matrix) = orthologous_genes
  # ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  ortholog_count_matrix[1:10,1:5]
  
  return(ortholog_count_matrix)
}

OrthologSeurat = function(speciesObject, orthology_key, common.genes = FALSE, mart_filepath = "../../Orthology/martMergeRefHuman.txt", reference_species = "Human"){
  
  species = GetSpecies(speciesObject)
  message("Working on ", species, "!")
  
  # Read in species data
  # if(species == "Mouse") {
  #   speciesObject = readRDS("../../Species_Reference/MouseACref_v4.rds")
  #   # speciesObject = subset(speciesObject, source == "Yan")
  # } else {
  #   speciesObject = readRDS(paste0(directory, species, "AC_v5.rds"))
  # }
  
  # Revert names to original species names
  if(!is.null(speciesObject@misc$orig.features)) {
    rownames(speciesObject@assays$RNA@counts) = speciesObject@misc$orig.features
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(species == "Rhabdomys") {
    # Using mouse as reference
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, "Mouse"))
  } else if(species == "Goldfish") {
    # Using zebrafish as reference
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, "Zebrafish"))
  } else if(!species %in% c(reference_species, "TreeShrew")) {
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, species))
  }
  
  # Get count matrix
  count_matrix = speciesObject@assays$RNA@counts
  
  # Now using same conversion function for consistency
  # # Convert some names to proper column name in orthology file
  # if(species == "Peromyscus") species = "Northern American deer mouse"
  # if(species == "Macaque") species = "Crab-eating macaque"
  # if(species == "Lizard") species = "Green anole"
  # if(species == "MouseLemur") species = "Mouse Lemur"
  # if(species == "Marmoset") species = "White-tufted-ear marmoset"
  # if(species == "Rhabdomys") species = "Mouse" # Rhabdomys assembly used mouse names
  # 
  # Subset to genes in orthology key
  # if(species == "Human") {
  #   subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key$`Gene name`),]
  # } else {
  #   subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key[[paste0(species," gene name")]]),]
  # }
  
  # Subset to genes in orthology key
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key$`Gene name`),]
  n_genes_found = nrow(subset_count_matrix)
  n_genes = length(unique(orthology_key$`Gene name`))
  message("Found ", n_genes_found, " genes out of ", n_genes, " total orthologs. Efficiency = ", n_genes_found/n_genes)
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  
  if(common.genes){
    ortholog_count_matrix = subset_count_matrix
  } else {
    # Now returning raw matrix with all genes, including genes that were not found in every species
    ortholog_count_matrix = subset_count_matrix[match(orthology_key$`Gene name`, rownames(subset_count_matrix)),]
    rownames(ortholog_count_matrix) = orthology_key$`Gene name`
    ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  }
  
  # Create new seurat object with metadata
  # speciesOrtho = CreateSeuratObject(subset_count_matrix)
  speciesOrtho = CreateSeuratObject(ortholog_count_matrix)
  speciesOrtho@meta.data = speciesObject@meta.data
  
  return(speciesOrtho)
}

# convertGenes <- function(genes, in_species, out_species = "Human"){
#   in_species_column=which(toupper(colnames(orthology_key))==paste0(species, " GENE NAME"))
#   orthologs = orthology_key[match(toupper(genes), toupper(orthology_key[,in_species_column])), paste0(out_species, " Gene Name")]
#   return(orthologs)
# }

ScoreCelltype = function(object, features){
  # Scale data if not already scaled for every feature
  # if(!all(rownames(object@assays$RNA@scale.data) == rownames(object))) object = ScaleData(object, features = rownames(object))
  features_in_object = intersect(features, rownames(object@assays$RNA@scale.data))
  message("Found ", length(features_in_object), " features in this object")
  raw_score = colSums(object@assays$RNA@scale.data[features_in_object,])
  normalized_score = raw_score/length(features_in_object)
  return(normalized_score)
}

AssignName = function(object, use = 'seurat_clusters', zero.index = FALSE){
  
  # Old function
  # speciesAC$type = paste0(sprintf('%02d', speciesAC$seurat_clusters), "_", speciesAC$lit_type) 
  # speciesAC$type[is.na(speciesAC$lit_type)] = sprintf('%02d', speciesAC$seurat_clusters[is.na(speciesAC$lit_type)])
  # speciesAC$type = factor(speciesAC$type, levels = 0:length(unique(speciesAC$type)))
  
  # object$type = paste0(str_pad(object$seurat_clusters, 2, side = "left"), "_", object$lit_type)
  # object$type[is.na(object$lit_type)] = str_pad(object$seurat_clusters[is.na(object$lit_type)], 2, side = "left")
  
  if(use == 'seurat_clusters' & zero.index) {
    object$clusters = as.numeric(as.character(object@meta.data[[use]]))
  } else if(use == 'seurat_clusters' & !zero.index) {
    object$clusters = as.numeric(as.character(object@meta.data[[use]]))+1
  } else {
    object$clusters = object@meta.data[[use]]
  }
  object$type = paste0(object$clusters, "_", object$lit_type)
  object$type[is.na(object$lit_type)] = as.character(object$clusters[is.na(object$lit_type)])
  meta = Metadata(object, feature.1 = "clusters", feature.2 = "lit_type")
  object$type = factor(object$type, levels = ifelse(is.na(meta$lit_type), as.character(meta$clusters), paste0(meta$clusters, "_", meta$lit_type)))
  object$type_no = object$clusters
  Idents(object) = "type"
  return(object)
}

nSpecies = function(object, group.by = "species"){
  return(length(unique(object@meta.data[[group.by]])))
}

nClusters = function(object, clusters = "seurat_clusters"){
  return(nrow(unique(object[[clusters]])))
}

GADPlot = function(df, Gaba.gene, Gly.gene = "SLC6A9", cutoff.Gly = NULL, cutoff.GABA = NULL){
  ggplot(df, 
         aes(x = eval(parse(text = Gaba.gene)), 
             y = eval(parse(text = Gly.gene)), 
             label = cluster)) + 
    geom_point()+
    xlab(Gaba.gene)+
    ylab(Gly.gene)+
    geom_hline(yintercept = cutoff.Gly, linetype = "dashed", color = "grey")+
    geom_vline(xintercept = cutoff.GABA, linetype = "dashed", color = "grey")+
    theme_dario() +
    scale_x_continuous(expand = expansion(mult = c(0,0)))+
    scale_y_continuous(expand = expansion(mult = c(0,0)))+
    geom_text_repel(min.segment.length = 0.05)
  
}

SubsampleSeurat = function(object, size){
  object <- object[, sample(colnames(object), size = size, replace=F)]
  return(object)
}

ExpressionScatter = function(data, x, y, xlab = NULL, ylab = NULL, title = NULL, 
                             genes.to.label = NULL, logX = FALSE, logY = FALSE, label.size = 4,
                             ...){
  
  if(inherits(data, "Seurat")) data = LogAvgExpr(data, assay = "RNA", ...)
  
  # Set all other genes to NA to avoid labeling
  data$gene = rownames(data)
  data$gene[!data$gene %in% c(genes.to.label)] = NA
  
  # Set labeled for coloring
  # data$labeled = ifelse(!is.na(data$gene), TRUE, FALSE)
  data$labeled = FALSE
  data$labeled[!is.na(data$gene)] = TRUE
  
  data$group = TRUE
  
  plt = ggplot(data, aes_(x= as.name(x), y=as.name(y), label = as.name("gene")))+#color = labeled, group = group))+
    ggrastr::rasterise(geom_point(data = subset(data, !labeled), color = "grey"), dpi = 300)+
    geom_point(data = subset(data, labeled), color = "black", shape = 21, fill = 'red')+
    # scale_color_manual(values = c("grey", "red"))+
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    {if(!is.null(genes.to.label)) geom_text_repel(max.overlaps = Inf, size = label.size)}+
    # labs(title = title,
    #       x = xlab,
    #      {if(!is.null(ylab)) y = ylab})+
    ggtitle(title)+
    {if(!is.null(xlab)) xlab(xlab)}+
    {if(!is.null(ylab)) ylab(ylab)}+
    stat_cor(method = "pearson", aes(label = ..r.label..))+
    {if(logX) scale_x_continuous(trans='log10')} +
    {if(logY) scale_y_continuous(trans='log10')} +
    theme_cowplot2()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
  return(plt)
}

#' @example
#' ScatterPlot(nTypes.df, 
#' x = "nCluster.residual", 
#' y = "nRgcCluster.residual", 
#' label = "species", 
#' ylab = "# of RGC clusters (residual)", 
#' xlab = "# of AC clusters (residual)", 
#' label.x.npc = 0.7, 
#' label.y.npc = 0.1)
ScatterPlot = function(data, x, y, xlab = NULL, ylab = NULL, title = NULL, 
                       labels = NULL, logX = FALSE, logY = FALSE, 
                       cor.position = c("top", "left"), 
                       lm = TRUE, r = FALSE, r_pval = TRUE, 
                       label_func = geom_text_repel, 
                       min.segment.length = 0.05, label.size = 3,
                       max.overlaps = 1, ...){
  
  if(!is.null(labels)) {
    data$names = data[[labels]]
  } else {
    data$names = ""
  }
  
  plt = ggplot(data, aes(x = eval(parse(text = x)), 
                         y = eval(parse(text = y)), 
                         label = names))+
    geom_point(fill = "grey", color = 'black', shape = 21)+
    labs(title = title, x = ifelse(is.null(xlab), x, xlab), y = ifelse(is.null(ylab), y, ylab))+
    {if(lm) geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")} +
    {if(r) stat_cor(method = "pearson", aes(label = ..r.label..), ...)}+
    {if(r_pval) stat_cor(method = "pearson", ...)}+
    {if(logX) scale_x_continuous(trans='log2')} +
    {if(logY) scale_y_continuous(trans='log2')} +
    {if(!is.null(labels)) label_func(max.overlaps = max.overlaps, min.segment.length = min.segment.length, size = label.size)}+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_dario()
  
  return(plt)
}

JSMatrix = function(matrix){
  rowSums = rowSums(matrix)
  colSums = colSums(matrix)
  
  js_matrix = matrix
  
  for(row in seq_along(rowSums)){
    for(column in seq_along(colSums)){
      intersection = matrix[row,column]
      union = rowSums[row] + colSums[column] - intersection
      js_matrix[row,column] = intersection/union
    }
  }
  
  return(js_matrix)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

JSHeatmap = function(matrix, title = NULL, stagger.threshold = 0.1, heatmap = TRUE, 
                     col.low = "white", col.high = "#584B9FFF", max.value = NULL,
                     legend_name = NULL, row.order = NULL, column.order = NULL, 
                     xlab = NULL, ylab = NULL, border.col = 'black',
                     split = FALSE){
  
  melted = reshape2::melt(matrix)
  colnames(melted) = c("row", "col", "jaccard")
  
  # Diagonalization
  row.max = apply(matrix, 1, which.max)
  match.df = data.frame(row = 1:length(row.max), 
                        col = as.integer(row.max))
  match.df$value = matrix[as.matrix(match.df)]
  match.df.sorted = match.df %>% arrange(-value)
  
  orig.rownames = rownames(matrix)
  orig.colnames = colnames(matrix)
  
  # If row.order is supplied, reorder the match.df.sorted
  if(!is.null(row.order)) match.df.sorted = match.df.sorted %>% arrange(factor(row, levels = match.df$row[match(row.order, orig.rownames)]))
  
  # Sorting
  melted.sorted = melted %>% arrange(factor(row, levels = unique(orig.rownames[match.df.sorted$row])), 
                                     factor(col, levels = unique(orig.colnames[match.df.sorted$col])))
  melted.sorted.nolow = subset(melted.sorted, jaccard > stagger.threshold)
  
  if(is.null(row.order)) {
    melted$row = factor(melted$row, levels = rev(unique(c(as.character(melted.sorted.nolow$row), setdiff(orig.rownames, melted.sorted.nolow$row)))))
  } else {
    melted$row = factor(melted$row, levels = rev(row.order))
  }
  
  if(is.null(column.order)) {
    melted$col = factor(melted$col, levels = (unique(c(as.character(melted.sorted.nolow$col), setdiff(orig.colnames, melted.sorted.nolow$col)))))
  } else {
    melted$col = factor(melted$col, levels = (column.order))
  }
  
  if(split) melted$species = str_split_fixed(melted$col, "_", 2)[,1]
  
  # Handle naming 
  if(is.null(legend_name)){
    if(heatmap) legend_name = 'Jaccard' else legend_name = 'Proportion'
  }
  
  if(is.null(max.value)){
    if(heatmap) max.value = 0.5 else max.value = 1
  }
  
  plt = ggplot(melted, aes(x = col, y = row))+
    {if(heatmap) geom_tile(aes(fill = jaccard), color=border.col) else geom_point(aes(colour = jaccard,  size=jaccard))}+
    {if(heatmap) scale_fill_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = col.high) else scale_color_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = col.high)} +
    # {{if(!heatmap) scale_size(range = c(1, max.size), limits = c(0, max.perc))+   
    theme_bw() +
    {if(heatmap) scale_y_discrete(expand = c(0,0))}+
    {if(heatmap) scale_x_discrete(expand = c(0,0))}+
    {if(split) facet_wrap(~species)}+
    {if(!heatmap) labs(size = legend_name)}+
    RotatedAxis() +
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(title)+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5))
  return(plt)
}

ConvertGeneSymbols = function(object, mart, species, make.unique = TRUE){
  
  if(inherits(mart, "character")){
    mart = readRDS(mart)
  }
  
  if(species == "Peromyscus") species = "Northern American deer mouse"
  if(species == "Macaque") species = "Crab-eating macaque"
  if(species == "Lizard") species = "Brown anole" #species = "Green anole"
  if(species == "MouseLemur") species = "Mouse Lemur"
  if(species == "Marmoset") species = "White-tufted-ear marmoset"
  
  if(species == "Lamprey"){
    lamprey_annotations = read.csv("../../../storage/data/GeneAnotation_Lamprey.xlsx - Sheet1.csv")
    lamprey_annotations$Symbol = gsub(" ", "", lamprey_annotations$Symbol)
    
    # Custom edits
    new_genes = rownames(object@assays$RNA@counts)
    # new_genes[new_genes == "LOC103091742"] = "TFAP2B_X2"
    # new_genes[new_genes == "MSTRG.3969"] = "TFAP2B_X3"
    
    # Convert from Gene first
    sorted_lamprey_annotations = lamprey_annotations[match(rownames(object@assays$RNA@counts), lamprey_annotations$Gene),]
    new_genes[new_genes %in% sorted_lamprey_annotations$Gene] = sorted_lamprey_annotations$Symbol[new_genes %in% sorted_lamprey_annotations$Gene]
    
    # Convert from Stringtie second
    sorted_lamprey_annotations = lamprey_annotations[match(rownames(object@assays$RNA@counts), lamprey_annotations$StringTie),]
    new_genes[new_genes %in% sorted_lamprey_annotations$StringTie] = sorted_lamprey_annotations$Symbol[new_genes %in% sorted_lamprey_annotations$StringTie]
    
    # Return blanks to original 
    new_genes[new_genes == ""] = rownames(object@assays$RNA@counts)[new_genes == ""]
    
    # Custom edits
    new_genes[new_genes == "SLC32A1/VIAAT1"] = "SLC32A1"
    new_genes[new_genes == "SLC6A9/GlyT1"] = "SLC6A9"
    new_genes[new_genes == "SLC6A1/GAT1"] = "SLC6A1"
    new_genes = make.unique(new_genes)
    
    # Print some examples of converted genes
    df = data.frame(before = rownames(object@assays$RNA@counts), after = new_genes)
    message("Converted ", nrow(df[df$before != df$after,]), " genes! Some examples below:")
    message(paste0(capture.output(head(df[df$before != df$after,], 200)), collapse = "\n"))
    
    # Add to object
    object@misc$orig.features = rownames(object@assays$RNA@counts)
    rownames(object@assays$RNA@counts) = new_genes
    rownames(object@assays$RNA@data) = new_genes 
  } else if(species == "TreeShrew"){
    TOGA_ts = read.table("../../Orthology/TOGA/TreeShrew.tsv", header = T)
    TOGA_ts$t_gene_name = str_split_fixed(TOGA_ts$t_transcript, "\\.", 2)[,2]
    TOGA_ts$q_gene_name = str_split_fixed(TOGA_ts$q_transcript, "\\.", 3)[,2]
    TOGA_ts_one2one = subset(TOGA_ts, orthology_class == "one2one")
    TOGA_ts_filt = unique(TOGA_ts_one2one[,c("q_gene","orthology_class","t_gene_name","q_gene_name")])
    new_genes = TOGA_ts_filt$t_gene_name[match(rownames(object), gsub("_", "-", TOGA_ts_filt$q_gene))]
    new_genes[is.na(new_genes)] = rownames(object)[is.na(new_genes)]
    rownames(object@assays$RNA@counts) = make.unique(new_genes)
    rownames(object@assays$RNA@data) = make.unique(new_genes)
  } else {
    # Pull genes to convert from object counts matrix
    before_genes = rownames(object@assays$RNA@counts)
    
    # Read in ensembl biomart orthology table
    # mart = as.data.frame(fread(mart_filepath))
    
    # Get one-to-one orthologs between species and mouse
    conversion_key = mart[mart[["Gene name"]] != "" & 
                            # mart[[paste0(species, " gene name")]] != "" & 
                            mart[[paste0(species, " homology type")]] == "ortholog_one2one",
                          c("Gene name", paste0(species, " gene stable ID"), paste0(species, " gene name"))] %>% unique
    
    # If the gene name is blank or not present in the object, try the ensembl id instead
    conversion_key[(conversion_key[[paste0(species, " gene name")]] == "" | !conversion_key[[paste0(species, " gene name")]] %in% before_genes), paste0(species, " gene name")] = 
      conversion_key[(conversion_key[[paste0(species, " gene name")]] == "" | !conversion_key[[paste0(species, " gene name")]] %in% before_genes), paste0(species, " gene stable ID")]
    
    # Sort genes by rownames of seurat object
    sorted_key = conversion_key[match(before_genes, conversion_key[[paste0(species, " gene name")]]),]
    
    # Convert seurat object names to reference (human or mouse) names if they are in the conversion key
    converted_genes = before_genes
    converted_genes[converted_genes %in% sorted_key[[paste0(species, " gene name")]] ] = toupper(sorted_key$`Gene name`[converted_genes %in% sorted_key[[paste0(species, " gene name")]] ])
    
    # In the rare case that we introduce duplicates, make unique to avoid gene merging and keep downstream analyses consistent
    if(make.unique) converted_genes = make.unique(converted_genes)
    
    # Change underscores to dashes
    converted_genes = gsub("_", "-", converted_genes)
    
    # Make conversion df
    before_after = data.frame(before = before_genes, after = converted_genes)
    
    # Print some examples of converted genes
    message("Converted ", nrow(before_after[converted_genes != before_genes,]), " genes! Some examples below:")
    message(paste0(capture.output(head(before_after[converted_genes != before_genes,], 200)), collapse = "\n"))
    
    # Save old symbols in the object
    object@misc$orig.features = before_genes
    
    # Change seurat object rownames
    rownames(object@assays$RNA@counts) <- converted_genes
    rownames(object@assays$RNA@data) <- converted_genes
    rownames(object@assays$RNA@scale.data) <- converted_genes[match(rownames(object@assays$RNA@scale.data), before_genes)]
  }
  return(object)
}

ProfileContamination = function(object, z.score.threshold = 3){
  object[["percent.rod"]] <- PercentageFeatureSet(object, features = intersect(Rod_markers, rownames(object)))
  object[["percent.cone"]] <- PercentageFeatureSet(object, features = intersect(Cone_markers, rownames(object)))
  object[["percent.hc"]] <- PercentageFeatureSet(object, features = intersect(HC_markers, rownames(object)))
  object[["percent.bc"]] <- PercentageFeatureSet(object, features = intersect(setdiff(BC_markers, "PRKCA"), rownames(object)))
  object[["percent.ac"]] <- PercentageFeatureSet(object, features = intersect(AC_markers, rownames(object)))
  object[["percent.rgc"]] <- PercentageFeatureSet(object, features = intersect(setdiff(RGC_markers, c("POU6F2", "RBFOX3")), rownames(object)))
  object[["percent.mg"]] <- PercentageFeatureSet(object, features = intersect(MG_markers, rownames(object)))
  
  # rod.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.rod")) > 2) - 1
  # cone.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.cone")) > 2) - 1
  # hc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.hc")) > 2) - 1
  # bc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.bc")) > 2) - 1
  # ac.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.ac")) > 2) - 1
  # rgc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.rgc")) > 2) - 1
  # mg.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.mg")) > 2) - 1
  
  cell_classes = c("rod", "cone", "hc", "bc", "ac", "rgc", "mg")
  
  contamin.clusters = lapply(cell_classes, function(class) {
    contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = paste0("percent.", class))) > z.score.threshold) - 1
    return(contamin)
  })
  
  names(contamin.clusters) = cell_classes
  object@misc$contamin.clusters = contamin.clusters
  return(object)
}

TrainXGBoost = function(train, test, train.clusters = "seurat_clusters", nfeatures = 2000){
  
  Idents(train) = train.clusters
  # Idents(test) = test.clusters
  
  DefaultAssay(train) = "RNA"
  DefaultAssay(test) = "RNA"
  
  train = FindVariableFeatures(train, selection.method = "vst", nfeatures = nfeatures)
  test = FindVariableFeatures(test, selection.method = "vst", nfeatures = nfeatures)
  common_HVGs <- intersect(VariableFeatures(train), VariableFeatures(test))
  message("Using ", length(common_HVGs), " common highly variable genes...\n")
  
  AC_model <- TrainModel(train, training_genes = common_HVGs, train_ident = train.clusters, do.scale = TRUE)
  
  return(AC_model)
}

LIBRARIES = c('paletteer', 'tidyverse', 'Seurat', 'ggplot2', 'reshape2', 
              'dplyr', 'xgboost', 'cowplot', 'pdfCluster', 'ggrepel', 
              'presto', 'scales', 'harmony', 'RColorBrewer', 'WGCNA', 
              'openxlsx', 'lisi', 'ggpubr', 'doParallel', 'circlize', 
              'viridis', 'ComplexHeatmap', 'Polychrome', 
              'HGNChelper', 'openxlsx', 'Matrix')

LoadLibraries = function(load.lisi = TRUE){

  res = lapply(LIBRARIES, require, character.only = TRUE)
  names(res) = LIBRARIES
  
  return(res)
}

SourceFiles = function(){
  source("../../utils/xgboost_train.R")
  source("../../utils/utilFxns.R")
  source("../../utils/xgboost_train.R")
  source("../../utils/plottingFxns.R")
  source("../../utils/dario_functions.R")
  source("../../utils/wrappers.R")
  # source("markers_metadata.R")
  source("../../utils/objects.R")
  source('../../utils/xgboost_train_DT.R')
  source('../../utils/SmartMatrix.R')
}

DoubletAnalysis = function(object, group.by = "seurat_clusters"){
  plot_grid(
    plot_grid(
      VlnPlot(object, features = "nFeature_RNA", group.by = "DF.classifications", pt.size = 0), 
      DimPlot(object, group.by = "DF.classifications") + NoLegend(), 
      DimPlot(object, group.by = group.by, label = T) + NoLegend(), 
      ncol = 3, nrow = 1, rel_widths = c(1, 2, 2)), 
    stackedBarGraph(object, feature.1 = "DF.classifications", feature.2 = group.by) + RotatedAxis(), 
    nrow = 2)
}

BrowseSeurat = function(object, batch = "animal"){
  if(is.null(object@meta.data[["percent.mt"]])) object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  plot = plot_grid(DimPlot(object, label = TRUE) + NoLegend(), 
                   ClusterBatchPlot(object, batch = batch, shuffle = TRUE) + NoLegend(), 
                   ClusterFeaturePlot(object, features = "nFeature_RNA") + NoLegend(),
                   VlnPlot(object, "nCount_RNA", pt.size = 0) + RotatedAxis() + NoLegend(), 
                   VlnPlot(object, "nFeature_RNA", pt.size = 0) + RotatedAxis() + NoLegend(), 
                   VlnPlot(object, "percent.mt", pt.size = 0) + RotatedAxis() + NoLegend(),
                   nrow = 2, ncol = 3, labels = LETTERS)
  return(plot)
}

BrowseSeuratSmall = function(object, batch = "orig.file"){
  plot = plot_grid(theme_umap(DimPlot(object, label = TRUE, group.by = "seurat_clusters", shuffle = TRUE)) + NoLegend(), 
                   theme_umap(ClusterBatchPlot(object, batch = batch, shuffle = TRUE)) + NoLegend(), 
                   theme_umap(ClusterFeaturePlot(object, features = "nFeature_RNA")) + NoLegend(),
                   nrow = 1, ncol = 3)
  return(plot)
}

ReprocessIntegrated = function(object, 
                               nPCs = 20, 
                               cluster_resolution = 1.5, 
                               nfeatures = 2000, 
                               selection.method = "vst", 
                               method = "seurat", 
                               k.param = 20, 
                               recompute.var.genes = FALSE, 
                               run.umap = TRUE, 
                               verbose = TRUE){
  if(method == 'seurat'){
    
    # Save old clusters
    object$old_seurat_clusters = object$seurat_clusters
    
    # Re-process using integrated assay
    DefaultAssay(object) = "integrated"
    if(recompute.var.genes) object = FindVariableFeatures(object, selection.method = selection.method, nfeatures = nfeatures)
    object <- ScaleData(object) %>% # vars.to.regress = "nCount_RNA" in ScaleData worsens batch correction somehow
      RunPCA() %>% 
      FindNeighbors(dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution, verbose = verbose)
    
    if(run.umap) object = RunUMAP(object, dims = 1:nPCs)
    
    DefaultAssay(object) = "RNA"
  } else if(method == "harmony"){
    
    # Save old clusters
    object$old_seurat_clusters = object$seurat_clusters
    
    object <- object %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
    if(run.umap) object = RunUMAP(object, reduction = "harmony", dims = 1:nPCs)
  } else {
    stop("method should be either 'seurat' or 'harmony'")
  }
  
  return(object)
}

ClusterBatchPlot <- function(object, batch = "orig.file", group.by = "seurat_clusters", ...){
  plot = DimPlot(object, group.by = batch, ...)
  plot$data$seurat_clusters <- as.factor(object@meta.data[[group.by]][match(rownames(plot$data), rownames(object@meta.data))])
  return(LabelClusters(plot, id = "seurat_clusters"))
}

ClusterFeaturePlot <- function(object, group.by = "seurat_clusters", ...){
  # args = list(...)
  plot = FeaturePlot(object, ...)
  plot$data$seurat_clusters <- as.factor(object@meta.data[[group.by]])
  return(LabelClusters(plot, id = "seurat_clusters"))
}

plotBatches2 = function(SeuratObject, batch = "orig.file", clusters = "seurat_clusters"){
  nGroups = length(unique(SeuratObject@meta.data[,batch]))
  p = plot_grid(
    ClusterBatchPlot(SeuratObject, batch = batch, group.by = clusters, shuffle = TRUE) + NoLegend(),
    plot_grid(
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch, split.by = batch) + NoLegend(), 
      stackedBarGraph(SeuratObject, feature.1 = batch, feature.2 = clusters) + RotatedAxis() + NoLegend(), 
      nrow = 2, 
      rel_heights = c(1,1.5)
    ),
    ncol = 2, 
    rel_widths = c(1, (1.5 + 0.1*nGroups))
  )
  return(p)
}

plotBatches = function(SeuratObject, batch = "orig.file", clusters = "seurat_clusters"){
  p = plot_grid(
    plot_grid(
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch) + NoLegend(),
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch, split.by = batch) + NoLegend(), 
      ncol = 2, 
      rel_widths = c(1.2, nrow(unique(SeuratObject[[batch]])))
    ),
    stackedBarGraph(SeuratObject, feature.1 = batch, feature.2 = clusters) + RotatedAxis(), 
    nrow = 2
  )
  return(p)
}

Harmonize = function(SeuratObject, batch = "orig.file", nPCs = 20, cluster_resolution = 0.5, k.param = 20, show.plots = TRUE, run.umap = TRUE, save.clusters = FALSE){
  
  # Save old clusters
  if(save.clusters) SeuratObject$old_seurat_clusters = SeuratObject$seurat_clusters
  
  DefaultAssay(SeuratObject) = "RNA"
  
  # Pre-processing
  SeuratObject <- SeuratObject %>% 
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(nPCs = nPCs, verbose = FALSE)
  
  harmony <- SeuratObject %>% RunHarmony(batch, plot_convergence = {if(show.plots) TRUE else FALSE})
  # Embeddings(mouse.harmony, 'harmony')[1:5, 1:5]
  
  # options(repr.plot.height = 5, repr.plot.width = 12)
  if(show.plots){
    print(plot_grid(
      DimPlot(SeuratObject, reduction = "pca", pt.size = .1, group.by = batch) + NoLegend(),
      DimPlot(harmony, reduction = "harmony", pt.size = .1, group.by = batch) + NoLegend()
    ))
  }
  
  # Downstream analysis
  if(!run.umap){
    harmony <- harmony %>%
      # RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
  } else {
    harmony <- harmony %>%
      RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
  }
  
  return(harmony)
}


plotClusterRes = function(SeuratObject, dims = 1:20, mc.cores = 10, from = 0.1, to = 2, by = 0.05){
  
  # Vary resolution parameter from 0.4 to 2.0
  resolutions = seq(from, to, by = by)
  library(parallel)
  nclusters = mclapply(resolutions, function(res){
    SeuratObject <- FindClusters(SeuratObject, resolution = res)
    return(length(unique(SeuratObject$seurat_clusters)))
  }, mc.cores = mc.cores) %>% unlist
  
  res.cluster.table = data.frame(resolution = resolutions,
                                 nclusters = nclusters)
  
  p = ggplot(res.cluster.table, aes(x = resolution, y = nclusters))+
    geom_line()+ 
    ylab("Number of clusters")+
    theme_bw()
  
  print(res.cluster.table)
  
  return(p)
}

proportionTable = function(SeuratObject, feature.1, feature.2){
  freq = table(SeuratObject@meta.data[[feature.1]], SeuratObject@meta.data[[feature.2]]) 
  prop = t(freq)/colSums(freq)
  return(prop)
}

#' Find doublets
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' @param SeuratObject Object of class Seurat
#' @param channels the identifier corresponding to the different 10X channels each cell was in
#' @param annotations cell type annotations for homotypic doublet estimation
#' @param suffix suffix in colnames(SeuratObject) that is often added automatically, e.g. ".1"
#' @param num.cores number of cores, default is 1
FindDoublets = function(SeuratObject, 
                        channels = "orig.file", 
                        annotations = "seurat_clusters", 
                        suffix = "", 
                        num.cores = 1, #length(unique(SeuratObject$orig.file)), 
                        classify.by = NULL, 
                        nPCs = 10){
  
  if("DF.classifications" %in% colnames(SeuratObject@meta.data)) SeuratObject$DF.classifications = NULL
  
  # Split object by sequencing channel since doublets can only arise between cells in the same channel
  obj.list <- SplitObject(SeuratObject, split.by = channels)
  
  # Run DoubletFinder on each sample separately
  library(DoubletFinder)
  obj.list = lapply(obj.list, function(object){
    
    # print(head(as.character(object@meta.data[[classify.by]])))
    
    # Pre-processing
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
    object <- ScaleData(object)
    object <- RunPCA(object, nPCs = nPCs)
    object <- RunUMAP(object, dims = 1:nPCs)
    
    # pK Identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(object, PCs = 1:nPCs, sct = FALSE, num.cores = 1)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    # BC maximization
    max_pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("Found max pK of: ", max_pK, "\n")
    
    # Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(object@meta.data[[annotations]]) ## ex: annotations <- object@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(object@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # Run DoubletFinder with varying classification stringencies
    # object <- doubletFinder_v3(object, PCs = 1:20, pN = 0.25, pK = max_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    if(!is.null(classify.by)){
      object <- doubletFinder_v3(object, 
                                 PCs = 1:nPCs, 
                                 pN = 0.25, 
                                 pK = max_pK, 
                                 nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE, 
                                 sct = FALSE, 
                                 annotations = as.character(object@meta.data[[classify.by]]))
    }
    else {
      object <- doubletFinder_v3(object, 
                                 PCs = 1:nPCs, 
                                 pN = 0.25, 
                                 pK = max_pK, 
                                 nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE, 
                                 sct = FALSE)
    }
    
    return(object)
  })
  
  # Extract doublet annotations
  extract_dfs = lapply(obj.list, function(object){
    columns = colnames(object@meta.data)[which(startsWith(colnames(object@meta.data), "DF."))]
    df = object@meta.data[columns]
    
    # Re-name
    if(!is.null(classify.by)){
      colnames(df) = c("DF.classifications", 
                       paste0(unlist(lapply(strsplit(columns[-1], "_"), function(x) x[[5]])), "_contribution"))
    } else {
      colnames(df) = "DF.classifications"
    }
    
    return(df)
  })
  names(extract_dfs) = NULL # to prevent rbind from appending list names to rownames
  extract_df = do.call(rbind, extract_dfs)
  
  # Save annotations for troubleshooting
  saveRDS(extract_df, "DF_annotations.rds")
  
  # rownames(extract_df) = unlist(lapply(strsplit(rownames(extract_df), "\\."), function(x) x[length(x)]))
  
  # Remove the prefix added from do.call(rbind)
  # prefix = paste0(unlist(lapply(strsplit(rownames(extract_df), "\\."), function(x) x[1])), ".")
  # prefix = paste0(unlist(lapply(names(obj.list), function(x) rep(x, nrow(obj.list[[x]])))), '.')
  # rownames(extract_df) = unlist(lapply(1:nrow(extract_df), function(i) sub(prefix[i], "", rownames(extract_df)[i])))
  
  # Add annotations to object
  if(!is.null(classify.by)){
    sorted_metadata = extract_df[match(colnames(SeuratObject), paste0(rownames(extract_df), suffix)),]
    SeuratObject = AddMetaData(SeuratObject, sorted_metadata)
  } else {
    SeuratObject$DF.classifications = extract_df[match(colnames(SeuratObject), paste0(rownames(extract_df), suffix)),]
  }
  
  return(SeuratObject)
}

#' Merge clusters fast!
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' 
quickMergeClusters = function(object, more.than.nDEGs = 6, log2FC.threshold = 1, p_val_adj.threshold = 0.001){
  tree <- object@tools$BuildClusterTree
  
  # Node heights
  library(ape)
  node.heights = as.integer(max(node.depth.edgelength(tree))-node.depth.edgelength(tree))
  names(node.heights) = 1:length(node.heights)
  order_nodes = sort(node.heights) #order(node.heights[-c(1:length(unique(YanAC$seurat_clusters)))])
  
  # Extract nodes below a given height (but not the tips), using 15% of max height as cutoff
  # nodes_to_check = names(order_nodes[order_nodes != 0 & order_nodes < height * max(node.heights)])
  nodes_to_check = head(names(order_nodes[order_nodes != 0]), 10)
  
  message("Checking nodes: ", paste0(nodes_to_check, collapse = ", "))
  
  node_nDEGs = lapply(nodes_to_check, function(node){
    node.markers <- wilcoxauc(object, ident.1 = 'clustertree', ident.2 = node)
    nDEGs(node.markers, log2FC.threshold = log2FC.threshold, p_val_adj.threshold = p_val_adj.threshold, presto = TRUE)
    # Same as doing: FindMarkers(YanAC, ident.1 = x, ident.2 = y)
  })
  
  message("number of DEGs: ", paste0(unlist(node_nDEGs), collapse = ", "))
  
  nodes_to_merge = nodes_to_check[which(node_nDEGs <= more.than.nDEGs)]
  
  message("Merging nodes: ", paste0(nodes_to_merge, collapse = ", "))
  
  for(node_to_merge in nodes_to_merge){
    merge_idents = sort(c(Seurat:::GetLeftDescendants(tree, node_to_merge) - 1, 
                          Seurat:::GetRightDescendants(tree, node_to_merge) - 1))
    
    # Set merged nodes to lower number of the two clusters
    object$seurat_clusters[WhichCells(object, idents = merge_idents)] = merge_idents[1] # paste0(merge_idents, collapse = "_")
  }
  
  # Make the clusters be in order
  old_clusters = sort(unique(object$seurat_clusters))
  new_clusters = 0:(length(old_clusters)-1)
  
  # Change each old cluster to its new cluster 
  for(i in 1:length(old_clusters)){
    object$seurat_clusters[object$seurat_clusters == old_clusters[i]] = new_clusters[i]
  }
  
  # Set the order and Idents
  object$seurat_clusters = factor(object$seurat_clusters, levels = new_clusters)
  Idents(object) = object$seurat_clusters
  
  return(object)
}

#' Merge clusters
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' 
mergeClusters = function(object, more.than.nDEGs = 6, log2FC.threshold = 1, p_val_adj.threshold = 0.001, num.nodes.to.check = 10){
  tree <- object@tools$BuildClusterTree
  
  # Node heights
  library(ape)
  node.heights = as.integer(max(node.depth.edgelength(tree))-node.depth.edgelength(tree))
  names(node.heights) = 1:length(node.heights)
  order_nodes = sort(node.heights) #order(node.heights[-c(1:length(unique(YanAC$seurat_clusters)))])
  
  # cat("Using ", height, " of ", max(node.heights), ": ", height * max(node.heights), "\n")
  # cat("Checking the ten most similar nodes...\n")
  
  # Extract nodes below a given height (but not the tips), using 15% of max height as cutoff
  # nodes_to_check = names(order_nodes[order_nodes != 0 & order_nodes < height * max(node.heights)])
  nodes_to_check = head(names(order_nodes[order_nodes != 0]), num.nodes.to.check)
  
  message("Checking nodes: ", paste0(nodes_to_check, collapse = ", "))
  
  node_nDEGs = lapply(nodes_to_check, function(node){
    node.markers <- FindMarkers(object, ident.1 = 'clustertree', ident.2 = node)
    nDEGs(node.markers, log2FC.threshold = log2FC.threshold, p_val_adj.threshold = p_val_adj.threshold)
    # Same as doing: FindMarkers(YanAC, ident.1 = x, ident.2 = y)
  })
  
  message("number of DEGs: ", paste0(unlist(node_nDEGs), collapse = ", "))
  
  nodes_to_merge = nodes_to_check[which(node_nDEGs <= more.than.nDEGs)]
  
  message("Merging nodes: ", paste0(nodes_to_merge, collapse = ", "))
  
  # order = levels(object$seurat_clusters)
  # object$seurat_clusters = as.character(object$seurat_clusters)
  
  for(node_to_merge in nodes_to_merge){
    merge_idents = sort(c(Seurat:::GetLeftDescendants(tree, node_to_merge) - 1, 
                          Seurat:::GetRightDescendants(tree, node_to_merge) - 1))
    
    # Set merged nodes to lower number of the two clusters
    object$seurat_clusters[WhichCells(object, idents = merge_idents)] = merge_idents[1] # paste0(merge_idents, collapse = "_")
  }
  
  # Make the clusters be in order
  old_clusters = sort(unique(object$seurat_clusters))
  new_clusters = 0:(length(old_clusters)-1)
  
  # Change each old cluster to its new cluster 
  for(i in 1:length(old_clusters)){
    object$seurat_clusters[object$seurat_clusters == old_clusters[i]] = new_clusters[i]
  }
  
  # for(element in 1:max){
  #   if(any(object$seurat_clusters == element)){
  #     next
  #   }
  #   object$seurat_clusters[object$seurat_clusters > element] = object$seurat_clusters[object$seurat_clusters > element] - 1
  #   cat(paste0("Re-numbering from ", element, "\n"))
  # }
  
  # Set the order and Idents
  object$seurat_clusters = factor(object$seurat_clusters, levels = new_clusters)
  Idents(object) = object$seurat_clusters
  
  return(object)
}

stackedBarGraph = function(SeuratObject, feature.1, feature.2, title = NULL){
  
  if(inherits(SeuratObject@meta.data[[feature.1]], "factor")) SeuratObject@meta.data[[feature.1]] = droplevels(SeuratObject@meta.data[[feature.1]])
  if(inherits(SeuratObject@meta.data[[feature.2]], "factor")) SeuratObject@meta.data[[feature.2]] = droplevels(SeuratObject@meta.data[[feature.2]])
  
  data <- table(SeuratObject@meta.data[[feature.1]], SeuratObject@meta.data[[feature.2]])
  melted = reshape2::melt(as.data.frame(data))
  # colnames(melted) = c("batch", "type", "composition")
  p = ggplot(melted, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    theme_minimal()+
    ggtitle(title)+
    xlab(feature.2)+
    labs(fill=feature.1)+
    theme_cowplot()+
    # scale_fill_manual(feature.1)+
    ylab('composition')+
    theme(axis.text.x = element_text(hjust = 1, angle = 45), plot.title = element_text(hjust = 0.5))
  return(p)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

nDEGs = function(table, up = TRUE, down = TRUE, 
                 log2FC.threshold = 0.25, 
                 p_val_adj.threshold = 0.05, 
                 presto = FALSE, value = FALSE){
  if(presto){
    if(up & !down){
      DEGs = (unique(subset(table, logFC > log2FC.threshold & padj < p_val_adj.threshold)$feature))
    } 
    else if(!up & down){
      DEGs = (unique(subset(table, logFC < -log2FC.threshold & padj < p_val_adj.threshold)$feature))
    } 
    else{
      DEGs = (unique(subset(table, abs(logFC) > log2FC.threshold & padj < p_val_adj.threshold)$feature))
    }
  }
  else {
    if(up & !down){
      DEGs = (unique(rownames(subset(table, avg_log2FC > log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    } 
    else if(!up & down){
      DEGs = (unique(rownames(subset(table, avg_log2FC < -log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    }
    else{
      DEGs = (unique(rownames(subset(table, abs(avg_log2FC) > log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    } 
  }
  
  if(value) return(DEGs) else return(length(DEGs))
}

#' Label transfer
#' 
#' Performs the Seurat label transfer procedure
#' 
#' @param reference reference with labels
#' @param query query to label
#' @param dims dimensions to use from PCA
#' @param threshold prediction score threshold for a label to be transferred
labelTransfer = function(reference, query, dims = 1:30, threshold = 0.7, ...){
  anchors <- FindTransferAnchors(reference = reference, 
                                 query = query, 
                                 dims = dims, 
                                 reference.reduction = "pca", 
                                 ...)
  predictions <- TransferData(anchorset = anchors, refdata = reference$transfer, dims = dims)
  query <- AddMetaData(query, metadata = predictions)
  
  # UMAP projection
  reference <- RunUMAP(reference, dims = dims, reduction = "pca", return.model = TRUE)
  query <- MapQuery(anchorset = anchors, reference = reference, 
                    query = query, refdata = list(transfer = "transfer"), 
                    reference.reduction = "pca", reduction.model = "umap")
  
  query$predicted.transfer[query$prediction.score.max < threshold] = "Unknown" # threshold based on figure 3D of Stuart et al. 2019
  return(query)
}

#' Clustered dotplot
#' 
#' Hierarchically clusters the columns of a DotPlot
#' 
#' @param dotplot a DotPlot that needs to be clustered
ClusteredDotplot = function(object, dendrogram = NULL, order = NULL, cluster.columns = TRUE, ...){
  
  if(!is.null(order)) {
    Idents(object) = order
    cluster.columns = FALSE
  }
  dotplot = DotPlot(object, ...) + coord_flip()
  
  df <- dotplot$data
  exp_mat<-df %>%
    dplyr::select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  row.names(exp_mat) <- exp_mat$features.plot  
  exp_mat <- exp_mat[,-1] %>% as.matrix()
  
  percent_mat<-df %>% 
    dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
    pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame() 
  row.names(percent_mat) <- percent_mat$features.plot  
  percent_mat <- percent_mat[,-1] %>% as.matrix()
  
  # col_fun = circlize::colorRamp2(c(-2, 0, 2), viridis(20)[c(1, 10, 20)])
  col_fun = circlize::colorRamp2(c(-2, 2), c("grey", "blue"))
  cell_fun = function(j, i, x, y, w, h, fill){
    # grid.rect(x = x, y = y, width = w, height = h,
    #           gp = gpar(col = NA, fill = NA))
    # grid.circle(x=x,y=y,r= percent_mat[i, j]/3 * min(unit.c(w, h)),
    #             gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
    # cell_fun = function(j, i, x, y, w, h, fill){
    #   grid.rect(x = x, y = y, width = w, height = h, 
    #             gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * 0.02,#unit.c(w),#min(unit.c(w, h)),
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
  
  if(is.null(dendrogram)){
    Heatmap(exp_mat,
            heatmap_legend_param=list(title="expression"),
            column_title = NULL, 
            col=col_fun,
            rect_gp = gpar(type = "none"),
            cell_fun = cell_fun,
            row_names_gp = gpar(fontsize = 10),
            # row_km = 4,
            cluster_rows = FALSE,
            cluster_columns = cluster.columns,
            cluster_row_slices = FALSE,
            border = "black", 
            use_raster = FALSE)
  } else {
    Heatmap(exp_mat,
            heatmap_legend_param=list(title="expression"),
            column_title = NULL, 
            col=col_fun,
            rect_gp = gpar(type = "none"),
            cell_fun = cell_fun,
            row_names_gp = gpar(fontsize = 10),
            # row_km = 4,
            cluster_columns = dendrogram,
            cluster_rows = FALSE,
            cluster_row_slices = FALSE,
            border = "black", 
            use_raster = FALSE)
  }
  
}

PrettyHistogram = function(vector, vline = NULL, logX = FALSE, logY = FALSE, title = NULL, xlab = NULL, bins = NULL, binwidth = NULL){
  # name = deparse(substitute(vector))
  df = data.frame(variable = vector)
  p = ggplot(df, aes(x=variable)) + 
    geom_histogram(color="black", fill="white", bins = bins, binwidth = binwidth)+ #aes(y=..density..))+
    xlab(xlab)+
    ylab("Frequency")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    {if(logX) scale_x_continuous(trans='log2')} +
    {if(logY) scale_y_continuous(trans='log2', expand = expansion(mult = c(0,0.1)))} +
    {if(!is.null(vline)) geom_vline(xintercept=vline, color="red3", linetype="dashed")}+
    theme_classic()+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

BINplot <- function(SeuratObject, x.limits=NULL, y.limits=NULL, title=NULL, 
                    feature.low = 200, feature.high = 7000, percent.mt.low = 0, percent.mt.high = 5){
  dat=data.frame(percent.mt=SeuratObject$percent.mt, 
                 count=SeuratObject$nCount_RNA, 
                 feature=SeuratObject$nFeature_RNA)
  plot=ggplot(dat, aes(percent.mt, feature))+
    geom_bin_2d(binwidth = c(0.1, 100))+
    scale_fill_gradient(low="grey", high="black", trans = "log")+
    geom_hline(yintercept=feature.low, color="red3", linetype="dashed")+
    geom_hline(yintercept=feature.high, color="red3", linetype="dashed")+
    geom_vline(xintercept=percent.mt.high, color="red3", linetype="dashed")+
    geom_vline(xintercept=percent.mt.low, color="red3", linetype="dashed")+
    {if(!is.null(x.limits)) scale_x_continuous(limits = x.limits)}+
    {if(!is.null(y.limits)) scale_y_continuous(limits = y.limits)}+
    {if(!is.null(title)) ggtitle(title)}+
    # stat_cor(method = "pearson", label.x = 20)+
    theme_classic()
  
  return(plot)
}



#' Convert genes 
#'
#' @param genes vector of genes
#' @param from the species from which the symbols originate
#' @param to the species to convert gene symbols to
#' 
#' @return Returns a vector of converted gene symbols
convertGenes <- function(genes, from, to){
  from_column=which(toupper(colnames(orthology_all_others))==paste0(toupper(from), ".GENE.NAME"))
  to_column=which(toupper(colnames(orthology_all_others))==paste0(toupper(to), ".GENE.NAME"))
  
  # If gene is not in orthology table, revert to original
  converted_genes = orthology_all_others[match(toupper(genes), toupper(orthology_all_others[,from_column])), to_column]
  converted_genes[is.na(converted_genes)] = genes[is.na(converted_genes)]
  
  return(converted_genes)
}

#' Annotated DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
AnnotatePlot = function(plt, annotation, 
                        object_x=NULL, x_annotation = NULL, color_genes_x = FALSE,
                        object_y=NULL, y_annotation = NULL, color_genes_y = FALSE, 
                        COLOR_FUN = ColorCode){
  
  y_axis_order = layer_scales(plt)$y$range$range
  x_axis_order = layer_scales(plt)$x$range$range
  
  annotate_y = (!is.null(y_annotation) | color_genes_y)
  annotate_x = (!is.null(x_annotation) | color_genes_x)
  
  if(annotate_y){
    if(!is.null(object_y)){
      y_idents = Idents(object_y) #@meta.data[,y_annotation]
      y_axis_annotation = object_y@meta.data[match(y_axis_order, y_idents),][,y_annotation]
      colors = COLOR_FUN(annotation, y_axis_annotation)
      # print(y_axis_annotation)
    } else if(color_genes_y){
      colors = Colors(annotation)[match(y_axis_order, Genes(annotation))]
    }
  } else if(annotate_x){
    if(!is.null(object_x)){
      x_idents = Idents(object_x) #@meta.data[,y_annotation]
      x_axis_annotation = object_x@meta.data[match(x_axis_order, x_idents),][,x_annotation]
      # print(x_axis_annotation)
      colors = COLOR_FUN(annotation, x_axis_annotation)
    }
    else if(color_genes_x){
      # x_axis_annotation = Annotation(annotation_object)[match(x_axis_order, Genes(annotation_object))]
      colors = Colors(annotation)[match(x_axis_order, Genes(annotation))]
    }
  } else {
    stop('Please add only an x or y axis annotation')
  }
  
  message("\nColors: \n", paste0(colors, collapse = ", "))
  
  plt2 = plt + 
    {if(annotate_y) annotate("rect",
                             ymin = seq(0.5, length(y_axis_order)-0.5, by = 1), ymax = seq(1.5, length(y_axis_order)+0.5, by = 1),
                             xmin = 0.5, xmax = length(x_axis_order)+0.5,
                             alpha = .4, fill = colors)}+
    {if(annotate_x) annotate("rect",
                             xmin = seq(0.5, length(x_axis_order)-0.5, by = 1), xmax = seq(1.5, length(x_axis_order)+0.5, by = 1),
                             ymin = 0.5, ymax = length(y_axis_order)+0.5,
                             alpha = .4, fill = colors)} +
    theme(panel.grid = element_line(color = rgb(235, 235, 235, 100, 
                                                maxColorValue = 255), 
                                    linewidth = 1, 
                                    linetype = 1))
  return(plt2)
}

#' Metadata from SeuratObject
#' 
Metadata = function(object, feature.1, feature.2, feature.3 = NULL, feature.4 = NULL){
  df = unique(object@meta.data[,c(feature.1, feature.2, feature.3, feature.4)]) %>% arrange(!!sym(feature.1))#arrange(eval(parse(text=feature.1)))
  rownames(df) = NULL
  return(df)
}

DEGsPerCluster = function(object, DE_table, group.by = "seurat_clusters"){
  DEGs = lapply(Clusters(object, group.by = group.by), function(cluster) {
    subset(DE_table, group == cluster & abs(logFC) > 1 & padj < 0.001)$feature
  })
  
  nDEG_data = data.frame(cluster = Clusters(object), 
                         nDEGs = unlist(lapply(DEGs, length)), 
                         size = as.numeric(table(object@meta.data[,group.by])))
  
  return(nDEG_data)
}

Clusters = function(object, group.by = "seurat_clusters"){
  clusters = sort(unique(object@meta.data[[group.by]]))
  return((clusters))
}

#' Metadata from SeuratObject
#' 
MeanMetadata = function(object, feature.1, feature.2){
  df = object@meta.data[,c(feature.1, feature.2)] %>% 
    arrange(eval(parse(text=feature.1)))
  rownames(df) = NULL
  
  clusters = Clusters(object, group.by = feature.1)
  means = lapply(clusters, function(cluster){
    mean(df[df[,feature.1] == cluster,feature.2])
  }) %>% unlist
  
  names(means) = Clusters(object, group.by = feature.1)
  return(means)
}

#' Ordered DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
OrderedDotPlot = function(object, 
                          annotation,
                          group.by = NULL, 
                          color.clusters.by = NULL, 
                          order = NULL,
                          ...){
  
  args = list(...)
  
  if(is.null(group.by)){
    labels = Idents(object)
  } else {
    labels = object[[group.by]]  
  }
  
  if(is.null(order)){
    order = CelltypeOrder(annotation)
  }
  
  # Arrange in proper order
  # object@meta.data[,color.clusters.by] = as.character(object@meta.data[,color.clusters.by])
  cell_class_df = Metadata(object, group.by, color.clusters.by) %>% 
    arrange(factor(eval(parse(text=color.clusters.by)), levels = order))
  
  # Set idents to order of color labels so that plot is in correct order
  Idents(object) = factor(object@meta.data[,group.by], levels = (cell_class_df[,group.by]))
  
  plot <- DotPlot(object, col.min = -1, col.max = 2, ...) + 
    scale_color_gradient(low = "lightgrey", high = "#584B9FFF", limits = c(-1,2))+
    scale_radius(limits = c(0,100), range = c(0, 6))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    guides(color = guide_colorbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed"))+
    theme(panel.grid = element_line(color = rgb(235, 235, 235, 100, 
                                                maxColorValue = 255), 
                                    linewidth = 1, 
                                    linetype = 1), 
          # Italicize gene symbols
          axis.text.x = element_text(face = "italic"))
  
  return(plot)
}

#' Annotated DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
AnnotatedDotPlot = function(object, 
                            annotation,
                            color.genes = FALSE,
                            ...){
  
  args = list(...)
  
  Idents(object) = args$group.by
  
  # Convert color.clusters.by to character since it causes issues when it's a factor
  object@meta.data[,args$color.clusters.by] = as.character(object@meta.data[,args$color.clusters.by])
  # if(is.factor(object[[args$group.by]])) object[[args$group.by]] = as.character(object[[args$group.by]]); message("converting factor to character\n")
  
  ordered.plot = OrderedDotPlot(object, annotation, ...)
  if(color.genes){
    annotated.plot = AnnotatePlot(ordered.plot, annotation, color_genes_x = TRUE) # y_annotation = args$color.clusters.by)
  } else{
    annotated.plot = AnnotatePlot(ordered.plot, annotation, object_y=object, y_annotation = args$color.clusters.by)
  }
  
  return(annotated.plot)
}








# Moratorium
major_class_colors = function(annotations){
  for(i in 1:nrow(major_color_code)){
    annotations[annotations == major_color_code$annotation[i]] = major_color_code$color[i]
  }
  return(annotations)       
}

peptide_colors = function(annotations){
  for(i in 1:nrow(peptide_color_code)){
    annotations[annotations == peptide_color_code$annotation[i]] = peptide_color_code$color[i]
  }
  return(annotations)       
}

#' Color code amacrines
#' 
#' Color codes a vector of amacrine cell annotations with a consistent color code
#'
#' @param annotations a vector of amacrine cell annotations
#'
#' @return Returns a color vector
amacrine_colors = function(annotations){
  for(i in 1:nrow(amacrine_color_code)){
    annotations[annotations == amacrine_color_code$annotation[i]] = amacrine_color_code$color[i]
  }
  return(annotations)       
}

#' Amacrine order
#' 
#' A function that returns the amacrine cell default order
#' 
amacrine_order = function(){
  return(amacrine_color_code$annotation)
}

peptide_order = function(){
  return(peptide_color_code$annotation)
}

major_class_order = function(){
  return(c("Rod", "Cone", "HC", "BC", "BP", "GabaAC", "GlyAC", "RGC", "MG", "Other"))
}

paramSweep_DT <- function(seu, PCs=1:10, sct = FALSE, num.cores=6) {
  require(Seurat); require(fields); require(parallel);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    # cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v3,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,
                        mc.cores=num.cores)
    # stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep_v3,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

