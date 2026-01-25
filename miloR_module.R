# miloR
register(MulticoreParam(workers = 80, progressbar = TRUE))
seu %<>% UpdateSeuratObject()
seu@active.assay = 'RNA'
seu %<>% JoinLayers()

seu.bak = seu
as.obj = CreateAssayObject(counts = seu[['RNA']]$counts)
seu[["RNA"]] = as.obj
# seu@active.assay = 'RNA'
library(miloR)
sce = as.SingleCellExperiment(seu)
reducedDim(sce, "UMAP") = Embeddings(seu, "umap")
reducedDim(sce, "PCA") = Embeddings(seu, "pca")
sce.milo = sce %>% Milo()
reducedDim(sce.milo, "UMAP") = reducedDim(sce,"UMAP")
reducedDim(sce.milo, "PCA") = reducedDim(sce,"PCA")
sce.milo %<>% buildGraph(k = 15, d = 30, reduced.dim = "PCA")
sce.milo %<>% makeNhoods(prop = 0.2, k = 15, d = 30, 
                         refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sce.milo)  # suggested 50-100 at approx.100%

# make sure Grouping parameters are correct
sce.milo %<>% countCells(meta.data = as.data.frame(colData(sce.milo)),
                         sample = "Sample")
traj_design = data.frame(colData(sce.milo))[,c("Sample", "Group")]
traj_design = distinct(traj_design)
rownames(traj_design) = traj_design$Sample
traj_design = traj_design[colnames(nhoodCounts(sce.milo)), , drop=FALSE]
traj_design$Group = factor(traj_design$Group, level = c("P", "T"))
traj_design
sce.milo %<>% calcNhoodDistance(d = 30, reduced.dim = "PCA")
da_results = testNhoods(sce.milo, design = ~ Group, design.df = traj_design,
                        fdr.weighting = "graph-overlap", norm.method = "TMM")

sce.milo %<>% buildNhoodGraph()
da_results = annotateNhoods(sce.milo, 
                            da_results, 
                            coldata_col = "type")

plotDAbeeswarm(da_results, group.by = "type",alpha = 0.1) +
  scale_color_gradient(low = "#7DB391", high = '#1B3955')
