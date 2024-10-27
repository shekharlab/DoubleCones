# Wrappers for functions

InstallMyPackages = function(libraries = LIBRARIES){
  
  lapply(libraries, install_package)
  
}

RunOU = function(df, expression){
  options(warn = 2)
  
  tree <- with(df,ouchtree(node,ancestor,Time/max(Time),spcode))
  
  features = colnames(expression)
  results = do.call(rbind, lapply(features, function(gene) {
    
    # Add gene to df
    df$gene = as.numeric(expression[,gene])
    message('working on gene ', gene)
    
    # h1 <- brown(tree_edges_full_final['calb1'],tree)
    h2 <- tryCatch(hansen(df['gene'],tree,(df['OU.1']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    # h3 <- hansen(tree_edges_full_final['calb1'],tree,tree_edges_full_final['OU.2'],sqrt.alpha=1,sigma=1)
    h4 <- tryCatch(hansen(df['gene'],tree,(df['OU.3']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    h5 <- tryCatch(hansen(df['gene'],tree,(df['OU.4']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    
    # Check for non-convergence
    if(is.null(h2) | is.null(h4) | is.null(h5)){
      return(data.frame(gene = gene, 
                        logLik_null = NA, 
                        logLik_primate = NA, 
                        logLik_rodent = NA,
                        lrt_primate = NA, 
                        lrt_rodent = NA,
                        p_value_primate = NA, 
                        p_value_rodent = NA))
    }
    
    # Get the log-likelihoods of both models
    logLik_null <- logLik(h2)
    logLik_primate <- logLik(h4)
    logLik_rodent <- logLik(h5)
    
    # Perform the Likelihood Ratio Test
    lrt_primate <- 2 * (logLik_primate - logLik_null)
    lrt_rodent <- 2 * (logLik_rodent - logLik_null)
    
    # Calculate the p-value using chi-square distribution
    p_value_primate <- pchisq(lrt_primate, df = summary(h4)$dof - summary(h2)$dof, lower.tail = FALSE)
    p_value_rodent <- pchisq(lrt_rodent, df = summary(h5)$dof - summary(h2)$dof, lower.tail = FALSE)
    
    return(data.frame(gene = gene, 
                      logLik_null = logLik_null, 
                      logLik_primate = logLik_primate, 
                      logLik_rodent = logLik_rodent,
                      lrt_primate = lrt_primate, 
                      lrt_rodent = lrt_rodent,
                      p_value_primate = p_value_primate, 
                      p_value_rodent = p_value_rodent))
  }))
  
  return(results)
}

GlaucomaDETest = function(this.object, group.by, ident.1, ident.2, filename, overwrite = FALSE, top.n = 2, 
                          cols = conditions_palette6, volcano.plot = TRUE, method = 'seurat', 
                          avg_log2FC_cutoff = 0.25, p_val_adj_cutoff = 0.05){
  
  filepath = dirname(filename)
  prefix = gsub('.xlsx', '', basename(filename))
  types.use = names(which(!apply(as.data.frame.matrix(table(this.object$type, this.object@meta.data[[group.by]])), 1, function(x) any(x == 0))))
  
  if(overwrite | !file.exists(filename)){
    if(method == 'seurat'){
      
      de_table = do.call(rbind, lapply(types.use, function(this.type){
          
          # tryCatch to avoid cases where there are no DEGs returned
          tryCatch({
            de = FindMarkersFast(subset(this.object, type == this.type), group.by = group.by, ident.1 = ident.1, ident.2 = ident.2, p_val_adj_cutoff = 1, avg_log2FC_cutoff = 0)
            de$type = factor(this.type, levels = ALLTYPES)
            
            if(volcano.plot){
              suppressWarnings(dir.create(paste0(filepath, '/volcano_plots/', prefix)))
              pdf(paste0(filepath, '/volcano_plots/', prefix, '/', this.type, '.pdf'))
              print(volcanoPlot(de, max_fdr = 1e-100, labels = TRUE, max.overlaps = 20, size = 3))
              dev.off()
            }
            
            return(de %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff))
          }, error = function(e) {
            return(NULL)
          })
      }))
  } else if(method == 'DESeq2'){
    library(DESeq2)
    
    pseudobulk <- AggregateExpression(this.object, assays = "RNA", slot = 'counts', return.seurat = TRUE, group.by = c(group.by, "type", 'sample_core'))
    pseudobulk[['condition']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,1]
    pseudobulk[['number']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,2]
    pseudobulk[['lit_type']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,3]
    pseudobulk[['sample_core']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,4]
    pseudobulk[['type']] = paste0(pseudobulk$number, '_', pseudobulk$lit_type)
    
    de_table = do.call(rbind, lapply(types.use, function(this.type){
      
      # tryCatch to avoid cases where there are no DEGs returned
      tryCatch({
        de = RunDESeq2(subset(pseudobulk, type == this.type), ident.1, ident.2)
        de$type = factor(this.type, levels = ALLTYPES)
        
        if(volcano.plot){
          suppressWarnings(dir.create(paste0(filepath, '/volcano_plots/', prefix)))
          pdf(paste0(filepath, '/volcano_plots/', prefix, '/', this.type, '.pdf'))
          print(volcanoPlot(de, fc_cutoff = 1, max_fdr = 1e-50, labels = TRUE, max.overlaps = 20, size = 3, max_fc = 6))
          dev.off()
        }
        
        return(de)
      }, error = function(e) {
        return(NULL)
      })
    }))
    
  } else {
    stop('method should be either seurat or DESeq2!')
  }
    de_table = subset(de_table, abs(avg_log2FC) > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)
    write.xlsx(de_table, file = filename) # saveRDS(de_table, filename)
  } 
  
  de_table = read.xlsx(filename) # readRDS(filename)
  de_table$type = factor(de_table$type, levels = ALLTYPES)
  
  de_up = as.data.frame(table(subset(de_table, avg_log2FC > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)$type)) %>% mutate(Direction = ident.1)
  de_down = as.data.frame(table(subset(de_table, avg_log2FC < avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)$type)) %>% mutate(Direction = ident.2)
  
  p1 = ggbarplot(rbind(de_up, de_down) %>% setNames(c('Type', 'Frequency', 'Direction')), 
                 x = 'Type', y = 'Frequency', fill = 'Direction', position = position_stack()) +
    scale_fill_manual(values = cols)+
    # NoLegend()+
    RotatedAxis()
  
  # Get top genes
  de_table %>% 
    group_by(type) %>%
    arrange(-abs(avg_log2FC)) %>% 
    slice_head(n = top.n) %>%
    ungroup() -> top.genes
  
  # p2 = VlnPlot(this.object,
  #         features = top.genes$gene,
  #         group.by = 'type',
  #         split.by = group.by,
  #         pt.size = 0,
  #         stack = TRUE,
  #         flip = TRUE,
  #         cols = conditions_palette6[match(c(ident.1, ident.2), names(conditions_palette6))])
  
  # top DEG
  top_gene = (de_table %>% arrange(p_val_adj, avg_log2FC))[1,]
  p2 = TitlePlot(
       VlnPlot2(this.object,
                features = top_gene$gene,
                split.by = 'sample_core',
                group.by = group.by,
                idents = top_gene$type,
                cols = conditions_palette6,
                stack = FALSE,
                combine= TRUE,
                pt.size = 0.1), 
       title = paste0(top_gene$gene, ' in ', top_gene$type))

  return(list(p1, p2))
  
}

SplitJaccardMatrix = function(object, group.by = 'species', ident.1 = 'seurat_clusters', ident.2 = 'type', nrow = 1, ncol = NULL, ...){
  
  # if(inherits(object, 'Seurat')) {
    objectList = SplitObject(object, split.by = group.by)
    # dfList = lapply(objectList, function(object) object@meta.data)
  # } else {
    
  # }
  
  # Order as original 
  if(is.factor(object@meta.data[[group.by]])) objectList = objectList[levels(object@meta.data[[group.by]])]
  
  heatmapList = lapply(seq_along(objectList), function(index) {
    JSHeatmap2(objectList[[index]]@meta.data[[ident.1]], objectList[[index]]@meta.data[[ident.2]], 
              title = names(objectList)[[index]], 
              ...)
  })
  names(heatmapList) = names(objectList)
  
  if(is.null(ncol)) ncol = length(objectList)
  
  # Plot
  ggarrange(plotlist = heatmapList, 
            ncol = ncol, 
            nrow = nrow, 
            common.legend = TRUE, 
            legend = "right", 
            align = 'h')
}

BinaryTree = function(genes = NULL, masks = NULL, return.masks = FALSE, ...){
  
  if(!is.null(genes)){
    full_masks = lapply(matrix.list[genes], BinarizeExpression, return.vector = FALSE)
    names(full_masks) = genes
    masks = lapply(full_masks, function(full_mask) {
      new = apply(full_mask, 2, mean)
      new[new > 0.5] = 1
      new[new <= 0.5] = 0
      return(new)
    })
    names(masks) = genes
    
    # Remove failures and TFs with no variance
    variance = sapply(masks, var)
    filtered.masks = masks[!is.na(masks) & variance != 0]
    filtered.full_masks = full_masks[!is.na(masks) & variance != 0]
    
    if(return.masks) return(list(filtered.full_masks, filtered.masks))
    
  } else if(!is.null(masks)){
    filtered.masks = masks
  } else{
    stop("One of genes or masks must be defined!")
  }
  
  nontf.expr.mask = do.call(rbind, filtered.masks)
  colnames(nontf.expr.mask) = colnames(matrix.list[[1]])
  nontf.order = rownames(nontf.expr.mask)
  nontf.expr.mask = nontf.expr.mask[nontf.order, ]
  colnames(nontf.expr.mask) = Metadata(bc_ortho, feature.1 = "seurat_clusters", feature.2 = "NOG")$NOG
  bc.orthotypes = colnames(nontf.expr.mask)
  
  # Heatmap(t(nontf.expr.mask), 
  #         rect_gp = gpar(col = "grey", lwd = 1), col = c("white", "black"), 
  #         show_row_dend = TRUE, show_column_dend = TRUE, 
  #         clustering_distance_rows = "manhattan", 
  #         clustering_method_rows = "average", 
  #         clustering_distance_columns = "manhattan", 
  #         clustering_method_columns = "average")
  
  nontf.trees = ConstructTrees(t(nontf.expr.mask), ...)
  
  return(nontf.trees)
  
}
  
ConstructTree = function(matrix, bootstrap = TRUE, seed = 42, root = NULL, nIter = 1000, mc.cores = 16){
  library(phangorn)
  library(TreeTools)
  
  # Set seed
  set.seed(seed)
  
  # Make phyDat object
  phydat = phyDat(matrix, type="USER", levels=c("0", "1"), compress = FALSE)
  
  # Distance methods
  treeUPGMA <- upgma(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  treeNJ <- NJ(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  
  # Maximum parsimony
  treeMP <- pratchet(phydat, trace = 0, minit=100)
  
  # Assign branch lengths
  treeMP  <- acctran(treeMP, phydat)
  treeMP  <- di2multi(treeMP, tol = 0) # Keep as binary tree
  if(inherits(treeMP, "multiPhylo")){
    treeMP <- unique(treeMP)
  }
  
  # Labeling internal nodes
  treeMP$node.label = (length(treeMP$tip.label)+1):((length(treeMP$tip.label)+treeMP$Nnode))
  
  if(!is.null(root)) {
    treeMP = root(phy = treeMP, node = root)
    # treeMP = midpoint(treeMP)
  } 
  
  return(list(phydat = phydat, UPGMA = treeUPGMA, NJ = treeNJ, MP = treeMP))
}

ConstructTrees = function(matrix, bootstrap = TRUE, seed = 42, root = NULL, nIter = 1000, mc.cores = 16){
  library(phangorn)
  library(TreeTools)
  
  # Set seed
  set.seed(seed)
  
  # Make phyDat object
  phydat = phyDat(matrix, type="USER", levels=c("0", "1"), compress = FALSE)
  
  # Distance methods
  treeUPGMA <- upgma(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  treeNJ <- NJ(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  
  # Maximum parsimony
  treeMP <- pratchet(phydat, trace = 0, minit=100)
  
  # Assign branch lengths
  treeMP  <- acctran(treeMP, phydat)
  treeMP  <- di2multi(treeMP, tol = 0) # Keep as binary tree
  if(inherits(treeMP, "multiPhylo")){
    treeMP <- unique(treeMP)
  }
  
  # Labeling internal nodes
  treeMP$node.label = (length(treeMP$tip.label)+1):((length(treeMP$tip.label)+treeMP$Nnode))
  
  if(!is.null(root)) {
    treeMP = root(phy = treeMP, node = root)
    # treeMP = midpoint(treeMP)
  } 
  
  # Bootstrap
  if(bootstrap){
    bs_upgma = as.multiPhylo(lapply(seq_len(nIter), function(iter) upgma(dist(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], method = "manhattan"))))
    treeUPGMA = plotBS(treeUPGMA, bs_upgma, main="UPGMA", type = "none")
    treeUPGMA$bs.support = treeUPGMA$node.label
    
    bs_nj = as.multiPhylo(lapply(seq_len(nIter), function(iter) NJ(dist(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], method = "manhattan"))))
    treeNJ = plotBS(treeNJ, bs_nj, main="NJ", type = "none")
    treeNJ$bs.support = treeNJ$node.label
    
    bs_mp = as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
      phydat = phyDat(matrix[,sample(1:ncol(matrix), replace = TRUE)], type="USER", levels=c("0", "1"), compress = FALSE)
      pratchet(phydat, trace = 0, minit=100)
    }, mc.cores = mc.cores))
    treeBS = plotBS(treeMP, bs_mp, main="MP", type = "none")
    treeMP$bs.support = treeBS$node.label
    # treeMP$bs.support[is.na(treeMP$bs.support)] = 0
  }
  
  return(list(phydat = phydat, UPGMA = treeUPGMA, NJ = treeNJ, MP = treeMP))
}

RunSoupx = function(object, feature.low = 800){
  
  # Find preliminary clusters for SoupX cleaning
  cells=ClusterSeurat(subset(object, nFeature_RNA >= feature.low))
  droplet_matrix = subset(object, nFeature_RNA < feature.low)@assays$RNA@counts
  
  # Make SoupChannel object
  sc = SoupChannel(droplet_matrix, # droplets
                   cells@assays$RNA@counts) # cells
  
  # Set clusters
  sc = setClusters(sc, cells$RNA_snn_res.0.5)
  
  # Automated method
  par(mfrow=c(1,3))
  sc = autoEstCont(sc)
  
  # Evidence of ambient RNA contamination in 1363, setting a higher contamination rate
  # sc.1363 = setContaminationFraction(sc.1363, 0.2)
  
  # save adjustment
  adjust = adjustCounts(sc)
  
  # Quantify corrected counts per cell
  cells$soupx.correction = colSums(cells@assays$RNA@counts-adjust)
  
  # Sanity check: verify that the appropriate fractions of counts were adjusted
  message(sum(cells$soupx.correction)/sum(cells@assays$RNA@counts))
  
  # Adjust counts
  cells@assays$RNA@counts=adjust
  
  # Combine
  cells = ClusterSeurat(cells)
  
  # Check SoupX correction
  umap.line = DimPlot(cells, reduction = "umap", group.by = "orig.ident")
  umap.soupx = FeaturePlot(cells, reduction = "umap", feature="soupx.correction")
  # df = data.frame(corrected.counts = cells$soupx.correction, orig.ident = cells$orig.ident)
  # violin = ggplot(df, aes(orig.ident, corrected.counts)) + geom_violin() + theme_bw()
  violin = VlnPlot(cells, group.by = 'orig.ident', features = 'soupx.correction')
  
  print(umap.line | umap.soupx | violin)
  
  return(cells)
}

OrthotypeHeatmaps2 = function(object){
  
  speciesACList = lapply(unique(object$species), function(currentSpecies) subset(object, species == currentSpecies))
  names(speciesACList) = unique(object$species)
  
  # Generate jaccard matrices
  jaccardList = lapply(seq_along(speciesACList), function(index) JSMatrix(t(table(speciesACList[[index]]$annotated, speciesACList[[index]]$seurat_clusters))))
  jaccardList = lapply(jaccardList, as.data.frame.matrix)
  
  # Get the best match in each species 
  bestMatch = lapply(jaccardList, function(matrix) apply(matrix, 1, max))
  bestMatchDf = as.data.frame(t(do.call(rbind, bestMatch)))
  colnames(bestMatchDf) = unique(object$species)
  bestMatchDf$OrthoType = rownames(bestMatchDf)
  
  # Melt and plot
  melted = reshape2::melt(bestMatchDf)
  
  ## Confusion matrices, sorted by OrthoType conservation
  OT_order = levels(reorder(melted$OrthoType, melted$value, FUN = median, decreasing = TRUE))
  heatmapList = lapply(seq_along(speciesACList), function(index) {
    JSHeatmap(JSMatrix(t(table(speciesACList[[index]]$annotated, speciesACList[[index]]$seurat_clusters))), 
              heatmap = TRUE, 
              title = unique(object$species)[[index]], 
              row.order = OT_order, 
              stagger.threshold = 0.25)
  })
  
  # Plot
  ggarrange(plotlist = heatmapList, 
            ncol = length(heatmapList), 
            nrow = 1, 
            common.legend = TRUE, 
            legend = "right")
}


OrthotypeAnalysis2 = function(reference, 
                              objectList, 
                              homologyList, 
                              group.by = "annotated",
                              types.use = NULL, 
                              downsample = 100, 
                              cluster_resolution = 0.4, 
                              sample.basis = NULL,
                              ...){
  
  # Subset homologyList to genes present in each object
  homologyList = lapply(seq_along(homologyList), function(index){
    homology = homologyList[[index]][,intersect(colnames(homologyList[[index]]), rownames(objectList[[index]]))]
    return(homology[rowSums(homology)>0,])
  })
  
  # Generate a basis (1:1 orthologs across all species)
  basis = Reduce(intersect, lapply(homologyList, function(x) rownames(x)))
  if(!is.null(sample.basis)) basis = sample(basis, size = sample.basis)
  message("Using a basis of ", length(basis), " genes!")
  
  # Subset to genes and types in ref
  if(length(types.use) > 0){
    seurat1 = DownsampleSeurat(subset(reference, annotated %in% types.use), group.by = group.by, size = downsample)
  } else {
    seurat1 = DownsampleSeurat(reference, group.by = group.by, size = downsample)
  }
  seurat1 = SubsetSeuratGenes(seurat1, basis)
  
  # Subset to genes and types in each species
  objectList = lapply(seq_along(objectList), function(index){
    
    # Downsample object
    if(length(types.use) > 0){
      seurat2 = DownsampleSeurat(subset(objectList[[index]], annotated %in% types.use), group.by = group.by, size = downsample)
    } else {
      seurat2 = DownsampleSeurat(objectList[[index]], group.by = group.by, size = downsample)
    }
    
    # Transform expression
    homology = homologyList[[index]][basis,]
    seurat2.converted = CreateSeuratObject(homology %*% seurat2@assays$RNA@counts[colnames(homologyList[[index]]),])
    
    # Transfer metadata
    seurat2.converted = TransferMetadata(seurat2, seurat2.converted)
    
    return(seurat2.converted)
  })
  
  # Merge
  ortho = merge(seurat1, objectList)
  
  # Without integration
  # ortho = ClusterSeurat(ortho)
  # print(DimPlot(ortho, group.by = c("species", "annotated")))
  
  # With integration
  ortho = ClusterSeurat(ortho, integrate.by = "species", cluster_resolution = cluster_resolution, ...)
  # print(DimPlot(ortho, group.by = "species") | DimPlot(ortho, label = TRUE) + NoLegend() | DimPlot(ortho, group.by = "annotated", label = TRUE) + NoLegend())
  # print(OrthotypeHeatmaps(ortho))
  
  return(ortho)
}

OrthotypeHeatmaps = function(object, stagger.threshold = 0.2){
  
  OrthoObjectList = SplitObject(object, split.by = "species")
  heatmapList = lapply(seq_along(OrthoObjectList), function(index) JSHeatmap(JSMatrix((table(OrthoObjectList[[index]]$seurat_clusters, 
                                                                                             OrthoObjectList[[index]]$annotated))), 
                                                                             heatmap = TRUE, 
                                                                             title = names(OrthoObjectList)[[index]], 
                                                                             stagger.threshold = stagger.threshold))
  
  # Plot
  ggarrange(plotlist = heatmapList, 
            ncol = length(OrthoObjectList), 
            nrow = 1, 
            common.legend = TRUE, 
            legend = "right", 
            align = "h")
}

OrthotypeAnalysis = function(objectList, orthology_key, downsample = 200, types.remove = NULL){
  
  # Generate the ortholog seurat object for each species
  speciesOrthoList = lapply(objectList, function(object) OrthologSeurat(object, 
                                                                        orthology_key = orthology_key, 
                                                                        common.genes = TRUE, 
                                                                        mart_filepath = "../../Orthology/martRefChicken.csv", 
                                                                        reference_species = "Chicken"))
  
  # Subset each matrix to common genes
  common.genes = Reduce(intersect, sapply(speciesOrthoList, function(x) rownames(x)))
  message("Found ", length(common.genes), " common genes")
  speciesOrthoList = lapply(speciesOrthoList, function(x) SubsetSeuratGenes(x, features = common.genes))
  
  # Downsample each species type to equal number of cells
  speciesSubsetList = lapply(speciesOrthoList, function(object) {
    Idents(object) = "annotated"
    object = subset(object, cells = WhichCells(object, downsample = downsample, seed = 12345))
    return(object)
  })
  
  # Merge
  OrthoObject <- merge(speciesSubsetList[[1]], y = speciesSubsetList[2:length(speciesSubsetList)])
  
  # Keep features in at least 3 cells
  # OrthoObject <- CreateSeuratObject(GetAssayData(OrthoObject), min.cells = 3)
  
  # Add species name
  # OrthoObject$species = rep(c("Chicken", "Lizard", "Opossum"), unlist(lapply(speciesSubsetList, ncol)))
  
  # Add celltype information
  OrthoObject$species_class = as.character(unlist(lapply(speciesSubsetList, function(speciesAC) speciesAC$cell_class)))
  OrthoObject$species_cluster = as.character(unlist(lapply(speciesSubsetList, function(speciesAC) speciesAC$annotated)))
  
  # Remove certain types
  if(!is.null(types.remove)) {
    OrthoObject = subset(OrthoObject, species_cluster %in% types.remove, invert = TRUE)
    message("Removing ", paste0(types.remove, collapse = ", "))
  }
  
  message("Using the following types: ")
  print(table(OrthoObject$species_class))
  print(table(OrthoObject$species_class, OrthoObject$species))
  
  # Save object
  # saveRDS(OrthoObject, output_file_v1)
  
  # Tabulate species
  message("# of cells from each species: ")
  print(t(t(table(OrthoObject[["species"]]))))
  
  # Integrate by species
  OrthoObject = ClusterSeurat(OrthoObject, integrate.by = "species", cluster_resolution = 0.5)
  
  return(OrthoObject)
}

ClusterEnrichmentComparison = function(object, cell.threshold = 1000){
  marker = ifelse(any(grep("NEUN", object$enrichment, ignore.case = T)), "RBFOX3", ifelse(grep("CD90", object$enrichment, ignore.case = T), "THY1", stop("Couldn't find marker")))
  
  # Save original clustering 
  object$orig.clusters = object$seurat_clusters
  objectList = SplitObject(object, split.by = "enrichment")
  
  # Include only enrichment batches with more than N cells
  objectList = objectList[sapply(objectList, ncol) > cell.threshold]
  
  # Downsample to same number of cells
  min_cells = min(sapply(objectList, function(object) length(Cells(object))))
  message(paste0("Subsetting down to ", min_cells, " cells!"))
  downsampledList = lapply(objectList, function(object) object[,sample(colnames(object), min_cells, replace = FALSE)])
  
  # Cluster each one separately
  downsampledList = lapply(downsampledList, function(obj) {
    if(length(unique(obj$animal)) > 1) {
      Harmonize(obj, batch = "animal", cluster_resolution = 1.5, run.umap = FALSE, show.plots = FALSE)
    } else {
      ClusterSeurat(obj, cluster_resolution = 1.5, do.umap = FALSE)
    }
  })
  
  # Order based on marker expression
  marker.order = rev(names(sort(as.matrix(AvgExpr(object, features = marker, assay = "RNA", group.by = "seurat_clusters"))[1,])))
  
  # Comparison to original clustering
  p.list = lapply(seq_along(downsampledList), function(index){
    object = downsampledList[[index]]
    name = names(downsampledList)[index]
    
    # How many original clusters were retrieved above given threshold
    stats = OverlapStatistics(table(object$orig.clusters, object$seurat_clusters))
    nRetrieved = length(unique(subset(stats, log.p >= 10 & overlap >= 30)$ident1))
    message("Retrieved ", nRetrieved, " clusters from ", name)
    
    JSHeatmap(JSMatrix((table(object$orig.clusters, object$seurat_clusters))), 
              heatmap = TRUE, 
              title = paste0(name, " subset"),
              row.order = marker.order,
              stagger.threshold = 0.25) + NoLegend()
  })
  
  
  object$ordered = factor(object$seurat_clusters, levels = rev(marker.order))
  
  return(plot_grid(plotlist = c(p.list,
                         list(VlnPlot(object, group.by = "ordered", features = marker, pt.size = 0) + coord_flip() + NoLegend() + NoAxes())), 
            ncol = length(p.list)+1, align = "h", axis = "bt",  rel_widths = c(rep(1, length(p.list)),0.4)))
  
}

ModifyRetinalAtlases = function(species){
  
  # Initial filepath
  initial_filepath = paste0("../../Species_Objects/", species, "_initial.rds")
  
  # Read in full object
  species_initial = readRDS(initial_filepath)
  
  # Print object to show number of cells and features
  message(paste0(capture.output(species_initial), collapse = "\n"))
  
  if(species %in% c("Macaque", "Marmoset", "Peromyscus", "Squirrel",
                    "Opossum", "Cow", "Sheep", "Lizard", "Mouse", 
                    "Pig", "Chicken", "Ferret", "Zebrafish", 
                    "MouseLemur", "Rat", "Mouse", "Lamprey")){
    species_initial = ConvertGeneSymbols(species_initial, "../../Orthology/martMergeRefHuman.txt", species) 
    
    # Check that original feature names were saved correctly and that no genes were lost
    stopifnot(length(rownames(species_initial@assays$RNA@counts)) == length(species_initial@misc$orig.features))
    
  } else if(species %in% c("Goldfish")){
    # goldfish_key = fread("../../Orthology/ZF_LA_SB_symbol_for_Seurat.txt")
    # new_genes = goldfish_key$V3[match(rownames(object), goldfish_key$V7)]
    # object@misc$orig.features = rownames(object@assays$RNA@counts)
    # rownames(object@assays$RNA@counts) = new_genes
    # rownames(object@assays$RNA@data) = new_genes
    
    # Goldfish has whole genome duplication so all symbols are duplicated
    species_initial = ConvertGeneSymbols(species_initial, "../../Orthology/martMergeRefHuman.txt", species, make.unique = FALSE) 
    metadata = species_initial@meta.data
    
    # Aggregate duplicated genes
    mydat = as.data.frame(as.matrix(species_initial@assays$RNA@counts))
    mydat$gene = rownames(species_initial)
    mydat.sum <- aggregate(. ~ gene, data = mydat, sum)
    
    # Check that it worked
    mydat[mydat$gene == "rbfox1",30:40]
    mydat.sum[mydat.sum$gene == "rbfox1",30:40]
    
    # Transfer metadata
    rownames(mydat.sum) = mydat.sum$gene
    species_initial2 <- CreateSeuratObject(mydat.sum[,colnames(mydat.sum) != "gene"])
    species_initial2 = TransferMetadata(from = species_initial, to = species_initial2)
    species_initial = species_initial2
    rm(species_initial2)
    
    # species_initial$orig.file = species_initial$orig.ident
    
    # Process for downstream analysis
    species_initial = ClusterSeurat(species_initial, cluster_resolution = 0.5)
    
    # Change meta features or will cause issues downstream
    # species_initial[["RNA"]]@meta.features <- data.frame(row.names = rownames(species_initial[["RNA"]]))
  }
  
  if(species %in% c("Macaque")) species_initial$seurat_clusters = species_initial$annotated
  Idents(species_initial) = "seurat_clusters"
  species_initial = UpperCase_genes(species_initial)
  # DotPlot(species_initial, features = rev(c("PAX6", "TFAP2A", "TFAP2B", "CHAT", "SLC5A7"))) + coord_flip() + RotatedAxis()
  
  # Recover mislabeled amacrines
  if(species == "Human") species_initial$cell_class[species_initial$seurat_clusters %in% c(5, 14, 24)] = "AC"
  if(species == "Marmoset") species_initial$cell_class[species_initial$seurat_clusters %in% c(27)] = "Rod"
  if(species == "Peromyscus") species_initial$cell_class[species_initial$seurat_clusters %in% c(14,48)] = "AC"
  if(species == "Cow") species_initial$cell_class[species_initial$seurat_clusters %in% c(8,9,20)] = "AC"
  if(species == "Opossum") species_initial$cell_class[species_initial$seurat_clusters %in% c(1,29,30)] = "AC"
  if(species == "Lizard") species_initial$cell_class[species_initial$seurat_clusters %in% c(18)] = "AC"
  if(species == "Zebrafish") species_initial$cell_class[species_initial$seurat_clusters %in% c(9,18,24)] = "AC"
  
  # Fails for chicken and zebrafish because seurat_clusters are not always in the same cell class
  if(species == "Chicken"){
    
    # Convert to major cell class
    species_initial$cell_class = as.character(species_initial$cell_class)
    species_initial$cell_class[species_initial$cell_class == "GabaAC" | species_initial$cell_class == "GlyAC"] = "AC"
    species_initial$cell_class[species_initial$cell_class == "BP"] = "BC"
    # species_initial$cell_class[!species_initial$cell_class %in% Annotation(major_annotation)] = "Other"
    species_initial$cell_class[species_initial$cell_class == "MicroG"] = "Other"
    
    # species_initial$cell_class = ifelse(startsWith(as.character(species_initial$annotated), "AC-"), "AC", "Other")
    DotPlot(species_initial, features = rev(Genes(major_annotation))) + coord_flip()
  } else if(species %in% c("Zebrafish", "Goldfish")){
    DotPlot(species_initial, features = rev(Genes(major_annotation))) + coord_flip()
  } else {
    
    # Convert to major cell class
    species_initial$cell_class = as.character(species_initial$cell_class)
    species_initial$cell_class[species_initial$cell_class == "GabaAC" | species_initial$cell_class == "GlyAC"] = "AC"
    species_initial$cell_class[species_initial$cell_class == "BP"] = "BC"
    # species_initial$cell_class[!species_initial$cell_class %in% Annotation(major_annotation)] = "Other"
    species_initial$cell_class[species_initial$cell_class == "MicroG"] = "Other"
    
    # Check if cones and rods are in same cluster
    if(max(apply(table(species_initial$cell_class, species_initial$seurat_clusters), 2, function(x) length(which(x != 0)))) > 1){
      species_initial$cell_class[species_initial$cell_class == "Rod" | species_initial$cell_class == "Cone"] = "PR"
    }
    
    species_initial$cell_class = factor(species_initial$cell_class, levels = unique(major_annotation@annotation))
    
    # Default palette for DimPlot
    # palette = c("cyan", "chartreuse2", "red", "gold", "magenta", "grey", "blue", "darkred")
    # if(!"Cone" %in% unique(species_initial$cell_class)) palette = c("cyan", "chartreuse2", "gold", "magenta", "grey", "blue", "darkred")
    
    # Plot
    print(plot_grid(
    {if(species != "Macaque") DimPlot(species_initial, group.by = "cell_class", cols = myPalette(7))}, # c("darkred", "red", "gold", "chartreuse2", "cyan", "blue", "magenta", "grey")
    AnnotatedDotPlot(species_initial, major_annotation_custom, group.by = "seurat_clusters", color.clusters.by = "cell_class",
                     features = rev(Genes(major_annotation_custom)), color.genes = FALSE) + coord_flip() + RotatedAxis(),
    rel_widths = c(1,2.5)))
  }
  
  return(species_initial)
}

RunSubsampling = function(object, nPermutations = 50, fraction.use = 0.8, shuffle.types = FALSE, nCores = 1, method, group.by){
  
  # Sample cells
  subsample.size = floor(length(Cells(object))*fraction.use)
  rand.cells.use = lapply(seq_len(nPermutations), function(x) sample(colnames(object), size = subsample.size, replace=F))
  
  # Sample parameters
  rand.nPCs = sample(c(15:30), nPermutations, replace = TRUE)
  rand.k = sample(c(16:24), nPermutations, replace = TRUE)
  rand.res = sample(seq(1, 2, by = 0.1), nPermutations, replace = TRUE)
  
  clusters = lapply(seq_len(nPermutations), function(iteration){
    
    message("Running iteration #", iteration, " with the following params: \n nPCs: ", rand.nPCs[[iteration]], "\n k.param: ", rand.k[[iteration]], "\n res: ", rand.res[[iteration]])
    
    sample = suppressMessages(SubsampleCluster(object, 
                                               cells.use = rand.cells.use[[iteration]],
                                               nPCs = rand.nPCs[[iteration]], 
                                               k.param = rand.k[[iteration]], 
                                               resolution = rand.res[[iteration]], 
                                               shuffle.types = shuffle.types, 
                                               method = method, 
                                               group.by = group.by))
    
    return(list(original = sample$old_seurat_clusters, permuted = sample$seurat_clusters))
  })#, mc.cores = nCores)
  
  # Name by permutation
  names(clusters) = seq_len(nPermutations)
  
  return(list(clusters = clusters, 
              rand.cells.use = rand.cells.use, 
              rand.nPCs = rand.nPCs, 
              rand.k = rand.k, 
              rand.res = rand.res))
}

SubsampleCluster = function(seurat, cells.use, nPCs, k.param, resolution, recompute.var.genes = FALSE, shuffle.types = FALSE, method = "seurat", group.by = "animal") {
  
  # Subsample
  message("Subsampling to ", length(cells.use), " cells!")
  subsample = seurat[, cells.use]
  
  # Re-process
  if(method == "harmony"){
    subsample = Harmonize(subsample, 
                          batch = group.by, 
                          nPCs = nPCs, 
                          k.param = k.param, 
                          cluster_resolution = resolution, 
                          show.plots = FALSE, 
                          run.umap = FALSE)
  } else if(method == "seurat"){
    subsample = ReprocessIntegrated(subsample, 
                                    nPCs = nPCs,
                                    k.param = k.param,
                                    cluster_resolution = resolution,
                                    recompute.var.genes = recompute.var.genes,
                                    run.umap = FALSE,
                                    verbose = FALSE,
                                    method = method)
  }
  
  return(subsample)
}

# Species-level amacrine cell analysis
render_report = function(species, 
                         initial = TRUE, 
                         batch_int = FALSE, 
                         integrate_by = "animal",
                         harmony = FALSE, 
                         contamination_threshold = 6, 
                         contamination = NULL, 
                         nFeature_threshold = -2, 
                         doublet_finder = FALSE, 
                         manual_annotation = NULL, 
                         sac_annotation = FALSE, 
                         de_expression = TRUE, 
                         save = TRUE, 
                         output_file = NULL){
  
  if(!is.null(output_file)) {
    output_file = output_file
  } else {
    output_file = paste0("html_reports/Report-", species, ".html")
  }
  
  rmarkdown::render("Species_AC_analysis_v11.Rmd", 
                    params = list(
                      species = species,
                      initial = initial,
                      batch_int = batch_int,
                      integrate_by = integrate_by,
                      harmony = harmony,
                      contamination_threshold = contamination_threshold,
                      contamination = contamination,
                      nFeature_threshold = nFeature_threshold,
                      doublet_finder = doublet_finder,
                      manual_annotation = manual_annotation,
                      sac_annotation = sac_annotation,
                      de_expression = de_expression,
                      save = save
                    ),
                    output_file = output_file
  )
}

ScTypeAnnotation = function(object, gs_list){
  library(HGNChelper)
  
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # DB file
  # db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  # tissue = "Eye" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = object[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  # In case Seurat is used, it is either seurat[["RNA"]]@scale.data (default), seurat[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or seurat[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(object@meta.data[object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Other"
  
  object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    object@meta.data$customclassif[object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  return(object)
}

AnnotatedUmap = function(object, annotation, group.by = "seurat_clusters", color.clusters.by = "cell_class", 
                         rel_widths = c(1,2.5), plot.dotplot = TRUE, color.genes = FALSE, title = NULL, coord_flip = TRUE, 
                         plot.proportions = FALSE, umap.legend = FALSE, pretty.umap = FALSE, ...){
  
  Idents(object) = group.by
  orig.object = object # AnnotatedDotplot function is sensitive to the factor levels of color.clusters.by
  
  object@meta.data[,color.clusters.by] = factor(factor(object@meta.data[,color.clusters.by], levels = unique(Annotation(annotation))))
  colors = GetAnnotationColors(levels(object@meta.data[,color.clusters.by]), annotation)
  message("Using the following colors: ", paste0(colors, collapse = ", "))
  
  if(pretty.umap){
    umap = PrettyUmap2(object, group.by = color.clusters.by, cols = colors, show.legend = ifelse(umap.legend, TRUE, FALSE), ...) 
  } else {
    umap = TitlePlot(theme_umap(ClusterBatchPlot(object, batch = color.clusters.by, cols = colors, shuffle = TRUE, group.by = group.by, ...)), title = title) + 
      {if(!umap.legend) NoLegend()}
  }
  
  if(plot.dotplot & plot.proportions){
    plt = plot_grid(
      # UMAP
      umap,
      
      # Dotplot
      AnnotatedDotPlot(orig.object, annotation, group.by = group.by, color.clusters.by = color.clusters.by,
                       features = rev(Genes(annotation)), color.genes = color.genes) + 
        RotatedAxis() + 
        NoLegend() +
        {if(coord_flip) coord_flip()} +
        theme(axis.title = element_blank()),
      
      # Barplot
      CelltypeProportionBarplot(object, 
                                show.all = FALSE) + 
        scale_fill_manual(values = pr_palette) + 
        theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
        coord_flip(),
      rel_widths = rel_widths,
      ncol = 3,
      align = "h", 
      axis = "bt"
    )
  }
  else if(plot.dotplot){
    plt = plot_grid(
      umap,
      AnnotatedDotPlot(orig.object, annotation, group.by = group.by, color.clusters.by = color.clusters.by,
                       features = rev(Genes(annotation)), color.genes = color.genes) + 
        RotatedAxis() + 
        {if(coord_flip) coord_flip()} +
        theme(axis.title = element_blank()),
      rel_widths = rel_widths, 
      align = "h", 
      axis = "bt"
    )
  } else {
    plt = umap
  }
  
  return(plt)
}

PlotAmacrineTypes = function(object, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", group.by = "seurat_clusters", annotation = NULL){
  
  # annotation = eval(parse(text = annotation))
  
  order = arrange(Metadata(object, feature.1 = group.by, feature.2 = "classification", feature.3 = "lit_type"), 
                  classification, 
                  factor(lit_type, levels = unique(Annotation(annotation))))[[group.by]]
  
  plt = AnnotatePlot(
          OrderedDotPlot(object, annotation, group.by = group.by, color.clusters.by = group.by, order = order, 
                         features = rev(c(Gaba.gene.1, Gaba.gene.2, "SLC6A9", setdiff(Genes(annotation), "SLC6A9")))), 
          annotation, object_y = object, y_annotation = "lit_type"
          ) + coord_flip() + theme(axis.title = element_blank())
  
  return(plt)
}

GabaGlyClassification2 = function(object, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", Gly.gene = "SLC6A9", 
                                  mus1 = c(0.01,0.5), mus2 = c(0.01,0.5), sigmas1 = c(0.01,0.1), sigmas2 = c(0.01,0.1), log.transform = FALSE){
  # Idents(object) = "seurat_clusters"
  # GabaGly.df = as.data.frame(t((AvgExpr(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data"))))
  # GabaGly.df = as.data.frame(apply(GabaGly.df, 2, rescale))
  # GabaGly.df$cluster = rownames(GabaGly.df)
  
  Idents(object) = "seurat_clusters"
  GabaGly.df = if(log.transform){
    t(log1p(AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA))
  } else {
    t((AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA))
  }
  GabaGly.df = as.data.frame(apply(GabaGly.df, 2, rescale))
  GabaGly.df$cluster = rownames(GabaGly.df)
  
  # Average of GAD1 and GAD2 scaled, normalized expression
  if(Gaba.gene.1 %in% colnames(GabaGly.df) & Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = paste0(Gaba.gene.1, "n", Gaba.gene.2)
    GabaGly.df[[Gaba.gene]] = rescale((GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]])/2)
    # GabaGly.df[[Gaba.gene]] = rescale(GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]]) # rescale(log2(rescale((GabaGly.df$GAD1 + GabaGly.df$GAD2), to = c(1,10))))
  } else if(Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.2
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.2]]
  } else if(Gaba.gene.1 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.1
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.1]]
  }
  
  # if("orig.file" %in% colnames(object@meta.data)) if(any(grepl("treeshrew", object$orig.file))) Gaba.gene = "SLC6A1"
  
  # Run EM for two Gaussian mixture model
  # GAD.model <- mixtools::normalmixEM(GabaGly.df[[Gaba.gene]], lambda=c(0.5,0.5), mu = mus1) # The starting values of sigma are important for the GAD model
  # Glyt1.model <- mixtools::normalmixEM(GabaGly.df[[Gly.gene]], lambda=c(0.5,0.5), mu = mus2)
  
  GAD.model <- mixtools::normalmixEM(GabaGly.df[[Gaba.gene]], lambda=c(0.5,0.5), mu=mus1, sigma=sigmas1) # The starting values of sigma are important for the GAD model
  Glyt1.model <- mixtools::normalmixEM(GabaGly.df[[Gly.gene]], lambda=c(0.5,0.5), mu=mus2, sigma=sigmas2)
  
  GABA.df = cbind(GabaGly.df[,c(Gaba.gene, Gly.gene, "cluster")], ComputePosteriors(GabaGly.df[[Gaba.gene]], GAD.model))
  Glyt1.df = cbind(GabaGly.df[,c(Gaba.gene, Gly.gene, "cluster")], ComputePosteriors(GabaGly.df[[Gly.gene]], Glyt1.model))         
  
  # print(plot_grid(
  #   MixtureHistogram(GABA.df, Gaba.gene, GAD.model, intersects = FindCutpoint(GAD.model), legend.name = "p_2_xi", xlab = paste0(Gaba.gene, " expression"), plot.points = TRUE), 
  #   MixtureHistogram(Glyt1.df, Gly.gene, Glyt1.model, intersects = FindCutpoint(Glyt1.model), legend.name = "p_2_xi", xlab = paste0(Gly.gene, " expression"), plot.points = TRUE), 
  #   ncol = 2
  # ))
  
  # print(cowplot::plot_grid({if(Gaba.gene.1 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.1)}, 
  #                          {if(Gaba.gene.2 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.2)}, 
  #                          GADPlot(GabaGly.df, Gaba.gene = Gaba.gene, Gly.gene = Gly.gene, cutoff.Gly = FindCutpoint(Glyt1.model), cutoff.GABA = FindCutpoint(GAD.model)) + 
  #                            xlab(Gaba.gene), 
  #                          ncol = 3))
  
  # Find cutpoints
  cutoff.GABA = FindCutpoint(GAD.model)
  cutoff.Gly = FindCutpoint(Glyt1.model)
  
  # Plot
  plt = plot_grid(
    MixtureHistogram(GABA.df, Gaba.gene, GAD.model, intersects = FindCutpoint(GAD.model), legend.name = "p_2_xi", xlab = paste0(Gaba.gene, " expression"), plot.points = TRUE) + 
      NoLegend() + rremove("xlab") + rremove("x.text"), NULL, NULL,
    NULL, NULL, NULL,
    GADPlot(GabaGly.df, Gaba.gene = Gaba.gene, Gly.gene = Gly.gene, cutoff.Gly = cutoff.Gly, cutoff.GABA = cutoff.GABA) + xlab(Gaba.gene), NULL, 
    MixtureHistogram(Glyt1.df, Gly.gene, Glyt1.model, intersects = FindCutpoint(Glyt1.model), legend.name = "p_2_xi", xlab = paste0(Gly.gene, " expression"), plot.points = TRUE) + 
      coord_flip() + NoLegend() + rremove("ylab") + rremove("y.text"), 
    ncol = 3, nrow = 3, rel_widths = c(2, 0, 1), rel_heights = c(1, 0, 2), align = "hv", axis = "bt")
  
  nGnGs = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                   eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  Dual.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gly.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                       eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gaba.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  
  message(length(nGnGs), " nGnGs types found: ", paste0(nGnGs, collapse = ", "))
  message(length(Dual.types), " dual types found: ", paste0(Dual.types, collapse = ", "))
  message(length(Gly.types), " glycinergic types found: ", paste0(Gly.types, collapse = ", "))
  message(length(Gaba.types), " GABAergic types found: ", paste0(Gaba.types, collapse = ", "))
  
  object$classification = NA
  object$classification[object$seurat_clusters %in% nGnGs] = "nGnG"
  object$classification[object$seurat_clusters %in% Dual.types] = "Both"
  object$classification[object$seurat_clusters %in% Gly.types] = "Gly"
  object$classification[object$seurat_clusters %in% Gaba.types] = "GABA"
  
  # Add posterior probabilities to seurat object
  # object$GABA_posterior = 
  # print(GABA.df)
  
  # Order based on group
  object$ac.order = factor(object$seurat_clusters, 
                           levels = arrange(Metadata(object, feature.1 = "seurat_clusters", feature.2 = "classification", feature.3 = "lit_type"), classification, lit_type)[["seurat_clusters"]])
  
  return(list(object = object, plt = plt))
}

GabaGlyClassification = function(object, group.by = 'seurat_clusters', cutoff.Gly = 0.15, cutoff.GABA = 0.15, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", Gly.gene = "SLC6A9", ac.order = FALSE){
  
  Idents(object) = group.by
  GabaGly.df = t(AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA)
  GabaGly.df = as.data.frame(apply(GabaGly.df, 2, rescale))
  GabaGly.df$cluster = rownames(GabaGly.df)
  
  # Sum of GAD1 and GAD2 average, scaled, normalized expression
  if(Gaba.gene.1 %in% colnames(GabaGly.df) & Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = paste0(Gaba.gene.1, "n", Gaba.gene.2)
    GabaGly.df[[Gaba.gene]] = rescale(GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]]) # rescale(log2(rescale((GabaGly.df$GAD1 + GabaGly.df$GAD2), to = c(1,10))))
  } else if(Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.2
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.2]]
  } else if(Gaba.gene.1 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.1
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.1]]
  }
  
  nGnGs = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                   eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  Dual.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gly.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                       eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gaba.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  
  message(length(nGnGs), " nGnGs types found: ", paste0(nGnGs, collapse = ", "))
  message(length(Dual.types), " dual types found: ", paste0(Dual.types, collapse = ", "))
  message(length(Gly.types), " glycinergic types found: ", paste0(Gly.types, collapse = ", "))
  message(length(Gaba.types), " GABAergic types found: ", paste0(Gaba.types, collapse = ", "))
  
  print(plot_grid({if(Gaba.gene.1 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.1, Gly.gene = Gly.gene)}, 
                  {if(Gaba.gene.2 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.2, Gly.gene = Gly.gene)}, 
                  GADPlot(GabaGly.df, Gaba.gene = Gaba.gene, Gly.gene = Gly.gene, cutoff.Gly = cutoff.Gly, cutoff.GABA = cutoff.GABA), 
                  ncol = 3))
  
  object$classification = NA
  object$classification[object@meta.data[[group.by]] %in% nGnGs] = "nGnG"
  object$classification[object@meta.data[[group.by]] %in% Dual.types] = "Both"
  object$classification[object@meta.data[[group.by]] %in% Gly.types] = "Gly"
  object$classification[object@meta.data[[group.by]] %in% Gaba.types] = "GABA"
  
  # Order based on group
  if(ac.order) object$ac.order = factor(object$seurat_clusters, levels = arrange(Metadata(object, feature.1 = group.by, feature.2 = "classification", feature.3 = "lit_type"), classification, lit_type)[[group.by]])
  
  return(object)
}
