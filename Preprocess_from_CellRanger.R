rm(list = ls());gc();setwd('~/r_repo/');set.seed(23112647)
source('~/r_repo/start_up_toolbox.r')

# Pre_QC -----------------------------------------------------------------------
sample = list.files('original_data/');sample
# log file
log = data.frame()

# Phenom data
sample = list.files('original_data/')
pd = fread('pd.csv')
pd = pd[match(sample,pd$Sample,),]

Group = pd$Group
Sample = pd$Sample

## Read in
sce.list = list()
seu.list = list()
for(i in 1:length(sample)){
  print(paste0('Creating Seurat Object: ', sample[i]))
  # as.Seurat
  seu.obj = CreateSeuratObject(counts = Read10X(paste0('original_data/',sample[i],'/raw_feature_bc_matrix')),
                               project = sample[i],
                               min.cells = 3,
                               min.features = 100)
  seu.obj$Group = Group[i]
  seu.obj$Sample = Sample[i]
  seu.obj$Barcodes = colnames(seu.obj)
  seu.obj = RenameCells(seu.obj, add.cell.id = sample[i])
  
  sce = as.SingleCellExperiment(seu.obj)
  colData(sce) = DataFrame(seu.obj@meta.data)
  sce$Sample = sample[i]
  log[i,'Sample'] = sample[i]
  
  ## ---- Remove empty droplets with DropletUtils ----
  sce.eD = emptyDrops(counts(sce),
                      lower = sort(colSums(counts(sce)))[2], 
                      BPPARAM = MulticoreParam(80))
  # filter 
  sce = sce[,which(sce.eD$FDR <= 0.001)] # cutoff FDR <= 0.01
  # 0.01 is cutoff from qiming's cell paper; 0.001 is default cutoff 
  
  temp = as.numeric(summary(sce.eD$FDR <= 0.01))
  log[i,'Empty Droplet Counts'] = temp[2]
  log[i,'Empty Droplet Percentage'] = round(100 * temp[2] / (temp[2]+temp[3]), 2)
  
  ## ---- Identify doublets with scDblFinder ----
  sce = scDblFinder(sce, 
                    cluster = TRUE,
                    verbose = F,
                    BPPARAM = MulticoreParam(80))
  # remove doublets
  temp = table(sce$scDblFinder.class)
  log[i,'Doublet Counts'] = temp[2]
  log[i,'Doublet Percent'] = round(100 * temp[2]/(temp[1] + temp[2]), 2)
  sce = sce[,sce$scDblFinder.class == 'singlet']
  keep.cells = colnames(sce)
  seu.obj = subset(seu.obj, cells = keep.cells)
  
  seu.list[[i]] = seu.obj
}

write.csv(log, file = 'LOG-PreProcess.csv')

seu = merge(x = seu.list[[1]], y = seu.list[-1])

# Filtering
seu$percent.mt = PercentageFeatureSet(seu, pattern = '^MT-')
seu$percent.rp = PercentageFeatureSet(seu, pattern = '^RP[SL]')
# VlnPlot(seu, c('nCount_RNA','nFeature_RNA','percent.mt','percent.rp'), ncol = 4)

seu %<>% 
  subset(subset = nFeature_RNA > 500
         & nFeature_RNA < 7500
         & percent.mt < 10
         & percent.rp < 50
         & nCount_RNA < 50000)

# Remove noise
all_genes = rownames(seu)
mt_genes = grep("^MT-", all_genes, value = TRUE)
rp_genes = grep("^RP[SL]", all_genes, value = TRUE)
hsp_genes = grep("^HSP", all_genes, value = TRUE)
ensg_genes = grep("^ENSG", all_genes, value = TRUE)
genes_to_remove = unique(c(mt_genes, rp_genes, hsp_genes, ensg_genes))
genes_to_keep = setdiff(all_genes, genes_to_remove)
seu = subset(seu, features = genes_to_keep)

seu %<>% 
  JoinLayers() %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  RunHarmony(group.by.vars = 'orig.ident', max.iter = 100) %>% 
  RunUMAP(reduction = 'harmony', dims = 1:30) %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:30)
seu %<>% FindClusters(algorithm = 4, resolution = 1.0)
