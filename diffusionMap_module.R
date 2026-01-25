library(destiny)
dims_use = 1:30
X = Embeddings(seu, "harmony")[, dims_use, drop = FALSE]
dm = DiffusionMap(X, k = 50, n_eigs = 30)
dc = eigenvectors(dm)
colnames(dc) = paste0("DC_", seq_len(ncol(dc)))
dc = dc[Cells(seu), , drop = FALSE]
seu[["diffmap"]] = CreateDimReducObject(embeddings = dc,
                                        key = "DC_",
                                        assay = DefaultAssay(seu))
dc2 = scale(dc[, 1:20])        
seu[["diffmap.scaled"]] = CreateDimReducObject(embeddings = dc2,
                                               key = "DCs_",
                                               assay = "RNA")
DefaultAssay(seu) = "RNA"
