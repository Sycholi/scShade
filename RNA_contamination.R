## decontX ----
seu %<>% JoinLayers()
library(decontX)
decontX_results = decontX(seu@assays$RNA@layers$counts)
seu$Contamination = decontX_results$contamination
FeaturePlot(seu, 'Contamination', raster = FALSE)
seu %<>% subset(Contamination <= 0.1)

## SoupX ----
seu = load10X('cellranger/outs/folder')
seu = autoEstCont(seu)
out = adjustCounts(seu)
