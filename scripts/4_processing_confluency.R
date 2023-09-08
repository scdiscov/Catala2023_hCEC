#_______________________________________________________________________________
# subsetting and analysing confluency dataset
#_______________________________________________________________________________

setwd("C:/Users/groen/OneDrive/SCD/scd-projects/MER-PC-4")
#--------------------- settings ------------------------------------------------

dir.create('output/confluency/',recursive=T, showWarnings = FALSE)
# install monocle3, seurat-wrappers
# devtools::install_github("satijalab/seurat-wrappers")


#--------------------- subset dataset ------------------------------------------
# import seurat object with all cells

scd <- readRDS('seurat_all.rds')


# select confluency samples
Seurat::Idents(scd) <- scd$Medium
scd <- subset(scd, idents = "Confluency")


#--------------------- normalize and process -----------------------------------

PCchoice <- 19 # for confluency dataset

scd <- Seurat::FindVariableFeatures(scd,selection.method = "vst", nfeatures = 2000)
scd <- Seurat::ScaleData(scd)
scd <- Seurat::RunPCA(scd, features = Seurat::VariableFeatures(object = scd))
# scd <- Seurat::RunUMAP(scd, dims = 1:PCchoice)
scd <- harmony::RunHarmony (scd, group.by.vars = c('Library','Donor'),
                            # plot_convergence = TRUE,
                            dims.use = 1:PCchoice)
scd <- Seurat::RunUMAP(scd, reduction="harmony", dims = 1:PCchoice, verbose = F)

umap <- Seurat::DimPlot(scd, group.by = 'Library')
pdf('output/confluency/umap.pdf',width = 7,height = 8); print(umap); dev.off()


#--------------------- clustering analysis  ------------------------------------
# identify clusters in the dataset

k.param = 100
resolution = 0.2
algorithm = 1

scd <- Seurat::FindNeighbors(scd, reduction = "harmony", dims = 1:PCchoice,k.param=k.param)
scd <- Seurat::FindClusters(scd,algorithm=algorithm, resolution = resolution)

clusters.umap <- Seurat::DimPlot(scd, reduction="umap", group.by = "seurat_clusters",repel=T)
pdf('output/confluency/umap_clusters.pdf',width = 7,height = 8); print(clusters.umap); dev.off()


# identify markers per cluster
markers <- Seurat::FindAllMarkers(scd,assay = "RNA", slot="data", only.pos = TRUE, min.pct = 0.33,
                          max.cells.per.ident = 1000, logfc.threshold = -Inf, test="wilcox", return.thresh = 0.01)
saveRDS(markers, file='output/confluency/DE_markers.rds')
saveRDS(scd, file='seurat_confluency.rds')


#--------------------- pseudotime analysis  ------------------------------------
# trajectory and pseudotime analysis using monocle 3

scd <- readRDS('seurat_confluency.rds') ## to remove
Seurat::DefaultAssay(scd)<-"RNA"

# filter out cluster 3
scd <- subset(scd, subset = seurat_clusters  == 3, invert = T)

# filter out lowly expressed genes
filter_genes <- function (object, min.value=1, min.cells = 0, genes = NULL) {
  # parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("filter_genes"))]
  # object <- Seurat:::SetCalcParams(object = object, calculation = "filter_genes", ... = parameters.to.store)
  genes.use <- rownames(object[["RNA"]]@data)
  # genes.use <- rownames(object[["RNA"]]@counts)
  
  if (!is.null(genes)) {
    genes.use <- intersect(genes.use, genes)
    object[["RNA"]]@counts <- object[["RNA"]]@counts[genes.use, ]
    return(object)
  } else if (min.cells > 0) {
    # num.cells <- Matrix::rowSums(object@data > min.value)
    num.cells <- Matrix::rowSums(object[["RNA"]]@counts > min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    object[["RNA"]]@counts <- object[["RNA"]]@counts[genes.use, ]
    object[["RNA"]]@data <- object[["RNA"]]@data[genes.use, ]
    return(object)
  } else {
    return(object)
  }
}
scd <- filter_genes(scd, min.value = 1,min.cells = ncol(scd)*0.05)

# convert to monocle compatible cell dataset
cds <- SeuratWrappers::as.cell_data_set(scd, graph="rna_snn",default.reduction = "umap")
cds <- monocle3::cluster_cells(cds)
seurat_clusters <- cds@colData[, 'seurat_clusters']
names(seurat_clusters) <- rownames(scd[[]])
cds@clusters$UMAP$clusters <- seurat_clusters
rm(seurat_clusters)

# calculate trajectory
cds <- monocle3::learn_graph(cds)
cds <- monocle3::order_cells(cds)

pseudo <- monocle3::plot_cells(cds, color_cells_by = "pseudotime") #label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE, cell_size = 2)
pdf('output/confluency/pseudotime.pdf',width = 8,height = 7); print(pseudo); dev.off()
rm(pseudo)

# transfer pseudotime information to seurat object to allow gene expression plotting
scd$Pseudotime <- monocle3::pseudotime(cds, reduction_method = "UMAP")[rownames(scd[[]])]


#--------------------- save object ---------------------------------------------

saveRDS(scd, file='output/seurat_confluency.rds')


#--------------------- integration with native CEC's  --------------------------
# construct Seurat object with data from 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186433
# filtered counts and metadata per cell are deposited


# load atlas seurat object
atlas <- readRDS("./data/atlas_Endothelium.rds") # or any name given to the reconstructed seurat object
# atlas$Library <- gsub("_.*","", rownames(atlas[[]]))
atlas$Library <- paste0("Atlas_",atlas$Sample)
atlas$Dataset <- rep('Atlas', nrow(atlas[[]]))
atlas$seurat_clusters <- paste0('Atlas_', atlas$seurat_clusters)

# match both metadata
scd$Library <- paste0("CEC_",scd$Library)
scd$Sample <- paste0(scd$Donor,'-',scd$Timepoint)
scd$seurat_clusters <- paste0('CEC_', scd$seurat_clusters)
scd$Dataset <- rep('CEC', nrow(scd[[]]))

to_keep <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'seurat_clusters','Library', 'Sample', 'Dataset')
scd@meta.data <- scd@meta.data[, to_keep]
rm(to_keep)

# merge datasets and process
data.list <- list('atlas' = atlas, 'scd' = scd) 
combined <- merge(data.list[[1]],data.list[c(-1)])
rm(data.list)

PCchoice = 20

combined <- Seurat::NormalizeData(combined, normalization.method = 'LogNormalize', scale.factor = 10000,verbose = TRUE)
combined <- Seurat::FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- Seurat::ScaleData(combined)
combined <- Seurat::RunPCA(combined, features = Seurat::VariableFeatures(object = combined),nfeatures.print = 0,verbose = F)
combined <- harmony::RunHarmony (combined, group.by.vars = c('Library'), 
                            # plot_convergence = TRUE, 
                            dims.use = 1:PCchoice)
combined <- Seurat::RunUMAP(combined, reduction="harmony", dims = 1:PCchoice,verbose = F)
Seurat::DimPlot(combined, group.by = 'Dataset')
Seurat::DimPlot(combined, group.by = 'seurat_clusters', label = T)


# predicted atlas cell clusters using Seurat's cell label transfer

# calculate the umap model for confl dataset
scd <- Seurat::RunUMAP(scd, dims = scd@commands$RunUMAP.RNA.harmony$dims, reduction = "harmony", return.model = TRUE)
# scd <- RunUMAP(scd, dims = scd@commands$RunUMAP.RNA.pca$dims, reduction = "pca", return.model = TRUE)

# calculate anchors 
data.anchors <- Seurat::FindTransferAnchors(reference = scd, query = atlas, dims = scd@commands$RunUMAP.RNA.harmony$dims, reference.reduction = "harmony")

# transfer cell type labels
atlas <- Seurat::TransferData(anchorset = data.anchors, reference = scd, query = atlas,
                    refdata = list(cluster = "seurat_clusters"))

predicted_clusters <- data.frame(cbind("Cell_nr" = table(atlas$predicted.cluster), "Fraction" = round(prop.table(table(atlas$predicted.cluster))*100,1)))
openxlsx::write.xlsx(predicted_clusters, file = 'output/confluency/Predicted_atlas_clusters.xlsx', rowNames = T, colNames = T, overwrite = T)
