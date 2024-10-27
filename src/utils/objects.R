myPalette = function(n_colors, darken = 0.2){
  # Amacrine palette
  return(colorspace::darken(as.character(paletteer_c("grDevices::Spectral", n_colors), darken)))
}

DcPalette = function(n_colors, colors = c("violet", "deepskyblue", "chartreuse3", "orange",  "orangered","brown", 'grey')){
  # Double cone palette
  return(colorRampPalette(colors = colors)(n_colors))
}

options(ggplot2.discrete.colour = lapply(c(4, seq(6,100, by = 4)), myPalette), 
        ggplot2.discrete.fill = lapply(c(4, seq(6,100, by = 4)), myPalette))

setClass("SeuratAnnotation", 
         slots=list(cluster = "vector", 
                    count = "vector",
                    mean.feature = "vector",
                    mean.mt = "vector",
                    pct.doublet = "vector", 
                    classification = "vector",
                    lit.type = "vector", 
                    neuropeptides = "list",
                    proportion = "data.frame"
                    ))

SeuratAnnotation = function(SeuratObject, 
                            group.by = "seurat_clusters", 
                            peptide.list = as.list(NA)){
  
  metadata = Metadata(SeuratObject, group.by, "lit_type", "classification")
  
  cluster = metadata$seurat_clusters
  
  classification = metadata$classification
  
  pct.doublet = 1-(proportionTable(SeuratObject, 
                                feature.1 = "DF.classifications", 
                                feature.2 = group.by)[,"Singlet"])
  
  lit.type = metadata$lit_type
  
  # Extract cluster counts
  count = as.vector(table(SeuratObject[[group.by, drop = TRUE]]))
  
  # Extract mean nFeature
  mean.feature = MeanMetadata(SeuratObject, group.by, "nFeature_RNA")
  mean.mt = MeanMetadata(SeuratObject, group.by, "percent.mt")
  
  # Proportions is simply a cross-tabulation with each enrichment group
  enriched.table = as.data.frame.matrix(table(SeuratObject[[group.by, drop = TRUE]], SeuratObject[["enrichment", drop = TRUE]]))
  
  # Return object
  new("SeuratAnnotation", 
      cluster = cluster, 
      count = count,
      mean.feature = mean.feature,
      mean.mt = mean.mt,
      pct.doublet = pct.doublet, 
      classification = classification,
      lit.type = lit.type, 
      neuropeptides = peptide.list, 
      proportion = enriched.table/rowSums(enriched.table))
}

setMethod("show", "SeuratAnnotation", function(object) {
  
  # Collapse neuropeptides into a comma-delimited list
  neuropeptides = lapply(object@neuropeptides, function(x) {
      paste0(x, collapse = ", ") 
    }) %>% unlist
  
  df = data.frame(cluster = object@cluster, 
                  count = object@count,
                   mean.feature = object@mean.feature,
                   mean.mt = object@mean.mt,
                   pct.doublet = object@pct.doublet, 
                   classification = object@classification,
                   proportion = object@proportion,
                   lit.type = object@lit.type, 
                   neuropeptides = neuropeptides, 
                   check.names = FALSE)
  
  rownames(df) = NULL
  
  print(df)
})

# https://stackoverflow.com/questions/19612839/set-method-initialize-s4-class-vs-using-function
as.data.frame.SeuratAnnotation = function(object) {
  
  # Collapse neuropeptides into a comma-delimited list
  neuropeptides = lapply(object@neuropeptides, function(x) {
    paste0(x, collapse = ", ") 
  }) %>% unlist
  
  df = data.frame(cluster = object@cluster, 
                  count = object@count,
                  mean.feature = object@mean.feature,
                  mean.mt = object@mean.mt,
                  pct.doublet = object@pct.doublet, 
                  classification = object@classification,
                  proportion = object@proportion,
                  lit.type = object@lit.type, 
                  neuropeptides = neuropeptides, 
                  check.names = FALSE)
  
  rownames(df) = NULL
  
  return(df)
}

setClass("CelltypeAnnotation", 
         slots=list(gene = "vector", 
                    annotation = "vector", 
                    color = "vector"))

CelltypeAnnotation = function(gene, annotation, color, to.upper = TRUE, include.other = FALSE){
  
  # Save order if annotation is a factor
  if(is.factor(annotation)){
    order = levels(annotation)
    annotation = as.character(annotation)
  } else {
    order = annotation
  }
  
  # Add fake gene for other
  if(include.other){
    gene$Other = 'NOGENE'
    annotation = c(annotation, 'Other')
    order = c(order, 'Other')
    color = c(color, 'grey')
  }
  
  # Remove duplicate genes
  duplicates = getDuplicates(unlist(gene))
  gene = lapply(gene, function(x) x[!x %in% duplicates])
  
  # Check for duplicated genes 
  stopifnot(length(unique(unlist(gene))) == length(unlist(gene)))
  
  # Convert gene symbols to upper-case for seurat
  if(to.upper) gene = lapply(gene, toupper)
  
  # Replicate annotations by number of genes for each
  annotation = unlist(lapply(1:length(gene), function(i) {
    rep(annotation[i], length(gene[[i]]))
  }))
  
  # Replicate colors by number of genes for each
  color = unlist(lapply(1:length(gene), function(i) {
    rep(color[i], length(gene[[i]]))
  }))
  
  # Sort in order
  initial = data.frame(gene = unlist(gene), 
                       annotation = annotation, 
                       color = color)
  sorted = initial %>% arrange(factor(annotation, levels = order))
  
  # Return CelltypeAnnotation object
  new("CelltypeAnnotation", 
      gene = sorted$gene, 
      annotation = sorted$annotation,
      color = sorted$color)
}

setMethod("show", "CelltypeAnnotation", function(object) {
  print(data.frame(gene = object@gene, 
             annotation = object@annotation, 
             color = object@color))
})

as.data.frame.CelltypeAnnotation = function(object){
  data.frame(gene = object@gene, 
             annotation = object@annotation, 
             color = object@color)
}

# setGeneric("ColorCode", 
#            function(x, ..., verbose = TRUE) standardGeneric("myGeneric"),
#            signature = "x"
# )

ColorCode = function(object, vector){
  
  # Anything that is not annotated goes to NA
  vector[!vector %in% Annotation(object)] = NA
  
  for(i in 1:length(object@annotation)){
    vector[vector == object@annotation[i]] = object@color[i]
  }
  
  return(vector)
}

CelltypeOrder = function(object) unique(object@annotation)

setGeneric("Genes", function(x, ...) standardGeneric("Genes"))
setMethod("Genes", "CelltypeAnnotation", function(x, annotation = NULL) {
  if(!is.null(annotation)) x@gene[x@annotation == annotation] else x@gene
  })

setGeneric("Annotation", function(x) standardGeneric("Annotation"))
setMethod("Annotation", "CelltypeAnnotation", function(x) x@annotation)

setGeneric("Colors", function(x) standardGeneric("Colors"))
setMethod("Colors", "CelltypeAnnotation", function(x) x@color)

# Amacrine subclass annotation
amacrine_annotation = CelltypeAnnotation(
                        gene = list(c("Chat", "Megf11", "Slc18a3", "Slc5a7", "Sox2"),
                                    c("NEFH", "Prkca", "Sdk1"), 
                                    c("NOS1"),
                                    c("Ddc", "Chl1", "Arhgdig"),
                                    c("TH"), 
                                    c("Slc17a7"), 
                                    c("Gjd2", "Prox1", "Dab1", "Nfia", "Dner", "Calb2"), 
                                    c("Satb2", "Ebf3"), 
                                    c("Slc17a8", "Sdk2")),
                          annotation = c("SAC",  
                                         "A17", 
                                         "nNOS", 
                                         "CA1",
                                         "CA2",
                                         "VG1", 
                                         "A2",  
                                         "SEG",  
                                         "VG3"), 
                        color = myPalette(9))
                          # color = c("red",
                          #                "orange", 
                          #                "gold", 
                          #                "chartreuse2", 
                          #                "darkgreen", 
                          #                "cyan", 
                          #                "skyblue", 
                          #                "darkblue", 
                          #                "magenta"))
                                         
amacrine_annotation2 = CelltypeAnnotation(
  gene = list(c("Chat", "Slc5a7", "Sox2"),
              c("NEFH", "Prkca", "Sdk1"), 
              c("NOS1"),
              c("Ddc", "TH"),
              c("TPBG"), 
              c("CARTPT"),
              c("SCGN"),
              c("SCG2"),
              c("PYGB"),
              c("CALB2", "npy", "cck", "calb1"),
              c("Tac1"),
              c("PDGFRA"),
              c("VIP"),
              c("Slc17a7"), 
              c("Gjd2", "Prox1", "Dab1", "Nfia", "Dner"), 
              c("SYT2"),
              c("Satb2", "Ebf3"), 
              c("Slc17a8", "Sdk2"), 
              c("PVALB"), 
              c("CRH"), 
              c("GBX2")),
  annotation = c("SAC",  
                 "A17", 
                 "nNOS", 
                 "CA1",
                 "CA2",
                 "CART",
                 "Spiny",
                 "Sec",
                 "Wiry",
                 "SemiL",
                 "SubP",
                 "T45",
                 "VIP",
                 "VG1", 
                 "A2",  
                 "A8",
                 "SEG",  
                 "VG3", 
                 "KT2", 
                 'CRH', 
                 'NNgaba'), 
  color = c("red",
                 "orange", 
                 "gold", 
                 "chartreuse2", 
                 "darkgreen", 
                 "grey30",
                 "grey80",
                 "grey30",
                 "grey80",
                 "grey30",
                 "grey80",
                 "grey30",
                 "darkred",
                 "cyan", 
                 "skyblue", 
                 "blue",
                 "darkblue", 
                 "magenta", 
                 "grey80", 
                  "grey80",
                  "grey30"))

amacrine_annotation_clean = CelltypeAnnotation(
  gene = list(c("CHAT", "SLC5A7"),
              c("PRKCA", "SDK1"), 
              c("NOS1"),
              c("DDC", "TH"),
              c("TPBG"), 
              c("PDGFRA"),
              c("VIP"),
              c("CRH"),
              c("GBX2"),
              c("NPY"),
              c("CARTPT", 'SLC35D3'),
              c('RXRG'),
              c('MAF'),
              
              c("GJD2"), 
              c("NFIB"),
              c("SATB2", 'EBF3'), 
              c("EBF2"), 
              c("SLC17A8", "SDK2"), 
              c("ROBO3"),
              c("TRHDE")),
  annotation = c("SAC",  
                 "A17", 
                 "nNOS", 
                 "CA1",
                 "CA2",
                 "PDGFRA",
                 "VIP",
                 "CRH",
                 "NNgaba",
                 "NPY",
                 "SLC35D3",
                 'RXRG',
                 'MAF',
                 
                 "A2",  
                 "A8",
                 "SEG", 
                 "NNgly",
                 "VG3", 
                 'ROBO3', 
                 'TRHDE'
                 ),
  color = c(head(myPalette(25), 13), tail(myPalette(25), 7))
)

# Major retinal class annotation
major_annotation = CelltypeAnnotation(
  gene = list(c("RBPMS", "SLC17A6", "NEFL", "NEFM"), #"POU6F2", "RBFOX3"), 
              c('VSX2', 'SCGN2', 'OTX2', 'ISL1', 'GRM6', 'GRIK1'), #"PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", "SLC6A9"), 
              c("ONECUT1", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3"), 
              c("SAG", "PDC", "RHO"), 
              # c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                       "BC",  # chartreuse2
                       "AC", # cyan
                       "HC", # gold
                       "Cone",# orange
                       "Rod",  # red
                       # "PR", # red
                       "MG"), # magenta
                 levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
                  "chartreuse2", 
                  "cyan", 
                  "gold", 
                  "red", 
                  "darkred", 
                  # "red",
                  "magenta"))

major_annotation_pr = CelltypeAnnotation(
  gene = list(c("RBPMS", "SLC17A6", "NEFL", "NEFM"), #"POU6F2", "RBFOX3"), 
              c('VSX2', 'SCGN2', 'OTX2', 'ISL1', 'GRM6', 'GRIK1'), #"PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", "SLC6A9"), 
              c("ONECUT1", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            "gold", 
            "red",
            "magenta"))
                  
major_annotation_custom = CelltypeAnnotation(
  gene = list(c("SAG", "PDC", "RHO"),
              c("PDE6H", "CRX", "ARR3"), 
              c("ONECUT1", "LHX1", "CALB1"), 
              c("VSX2", "CABP5", "GRIK1"), 
              c("TFAP2A", "GAD2", "SLC6A9"), 
              c("RBPMS", "SLC17A6", "NEFL"),  
              c("SLC1A3", "RLBP1",  "APOE"), 
              c("Other")), 
  annotation = factor(c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG", "Other"),
                      levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG", "Other")),  
  color = c(myPalette(7), "grey"))

Rod_markers = Genes(major_annotation, annotation = "Rod")
Cone_markers = Genes(major_annotation, annotation = "Cone")
HC_markers = Genes(major_annotation, annotation = "HC")
BC_markers = Genes(major_annotation, annotation = "BC")
AC_markers = Genes(major_annotation, annotation = "AC")
RGC_markers = Genes(major_annotation, annotation = "RGC")
MG_markers = Genes(major_annotation, annotation = "MG")

phylogenetic_order = c("Human","Macaque","Marmoset","MouseLemur","TreeShrew",
                       "Mouse","Rhabdomys", "Rat","Peromyscus","Squirrel",
                       "Ferret","Pig","Cow","Sheep","Opossum","Chicken",
                       "Lizard", "Zebrafish", "Goldfish", "Lamprey")

phylogenetic_order2 = c("Human","Macaque","Marmoset","MouseLemur","TreeShrew",
                       "Mouse","Rhabdomys", "Rat","Peromyscus","Squirrel",
                       "Cow","Sheep","Pig","Ferret","Opossum","Chicken",
                       "Lizard", "Zebrafish", "Goldfish", "Lamprey")

SPECIESFILES = c(Human = "Human", Macaque = "Macaque", Marmoset = "Marmoset", 
                 MouseLemur = "MouseLemur", TreeShrew = "TreeShrew", Mouse = "Mouse",
                 Rhabdomys = "Rhabdomys", Rat = "Rat", Peromyscus = "Peromyscus", 
                 Squirrel = "Squirrel", Cow = "Cow", Sheep = "Sheep", 
                 Pig = "Pig", Ferret = "Ferret", Opossum = "Opossum", 
                 Chicken = "Chicken_reclustered", Lizard = "Lizard_ncbi", 
                 Zebrafish = "Zebrafish", Goldfish = "Goldfish", Lamprey = "Lamprey")

# Photoreceptor type annotation
pr_annotation = CelltypeAnnotation(
  gene = list(c("full"), 
              c("opn1mw4/opn1lw1"),
              c("ZEB2", "zeb2a", "zeb2b"),
              c("opn1sw1", "OPN1SW"), 
              c("opn1sw2", "OPN2SW", "LOC132767847"),
              c("opn1mw1", "opn1mw2", "opn1mw3", "opn1mw4", "OPN1MSW", "RH2", "OPN1MW", "LOC132773706"),
              c("opn1lw1", "opn1lw2", "OPN1LW", "LOC132767849", "ENSSTOG00000024701"),
              c("THRB", "thrb"),
              # c("principal"),
              c("MYLK", "mylk", "STBD1", "stbd1"),
              c("RHO", "rhol")),
  annotation = factor(c("full",
                        "opn1mw4/\nopn1lw1",
                        "singleCone",
                        "UV", 
                        "blue", 
                        "green", 
                        "red", 
                        "principal",
                        # 'principle',
                        "accessory",
                        'rod'), 
                      levels = rev(c("rod", "singleCone", "UV", "blue", "green", "red", "opn1mw4/\nopn1lw1", "principal", "accessory", "full"))),  
  color = c("tan", "gold", "white", "violet", "deepskyblue", "chartreuse3", "orangered", "brown", "orange", "lightgrey"), 
  to.upper = FALSE)

# Species image files
image.files = list.files("../../figures/animals/black", full.names = TRUE) %>% setNames(c('Chicken', 'Human', 'Lizard', 'Opossum', 'Rat', 'Zebrafish'))

major_annotation_pr_ze = CelltypeAnnotation(
  gene = list(c('isl2b', 'rbpms2b', 'robo2'), 
              c('vsx1', 'vsx2', 'otx2', 'isl1', "prkcab", "prkcaa", 'grik1a', 'grik1b'),
              c('pax6a', 'tfap2a', 'tfap2b', 'tfap2c', 'gad2', 'slc6a9'), 
              c("onecut2"), 
              c('rhol', 'opn1lw1','opn1mw1', 'opn1sw1','opn1sw2', 'arr3a', 'pde6h'),
              c("slc1a3b", "rlbp1b",  "apoe", "apoba")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            "gold", 
            "red", 
            "magenta"), to.upper = FALSE)

# pr_annotation = CelltypeAnnotation(
#   gene = list(c("PTPRD", "OPN1SW", "NRXN3", "RS1", "GNGT2"), 
#               c("THRB", "SOX5", "ACSL1", "OPN1MW", "OPN1LW")
#               c("RHO", "ESPRRB", "GNAT1", "PDE6B", "PDE6G", "PDC")),
#   annotation = factor(c("UV", 
#                         "Blue", 
#                         "Green", 
#                         ), # magenta
#                       levels = c("Rod", "ML_cone", "S_cone")),  
#   color = c("violet", "blue", "green", "red", "grey"))

# Consistent colors across plots; don't use globally just yet
# options(ggplot2.discrete.colour= list(dtColors()))


