options(stringsAsFactors=FALSE)

library(Seurat)
library(dplyr)
library(Matrix)
library(SCopeLoomR)
library(IRdisplay)
library(plyr)
library(patchwork)
library(ggplot2)

seuratObj <- readRDS("/lustre1/project/stg_00135/Cecilia/Data/ProjectedSatija/projected.both.rds")

seuratObj1 <- subset(seuratObj, subset = subsample == "d. Pre-Treatment P123RS")
table(seuratObj1$subsample)
loom <- build_loom(file.name = "redo.pre.p116RS.loom", 
                        dgem = as.matrix(seuratObj1[["RNA"]]@counts),
                        title = "STAT1 GoF Dataset",
                        genome="Human",
                        default.embedding = seuratObj1@reductions$ref.umap@cell.embeddings,
                        default.embedding.name = "ref.umap")
loom <- open_loom("/lustre1/project/stg_00135/Cecilia/Data/SCENIC/pre-p123rs_pySCENIC.loom", mode = "r+")


add_seurat_clustering(loom, 
                      seurat = seuratObj1,
                      seurat.assay = "RNA",
                      default.clustering.resolution = 0.2,
                      seurat.clustering.prefix = "RNA_snn_res.")



for(keyName in c('subsample', 'RNA_MERGED_clusters','predicted.celltype.l1', 'predicted.celltype.l2'))
{
    message('Adding ', keyName, '...')
    add_col_attr(loom, key = keyName,
                value = as.vector(seuratObj1[[keyName]][,1]),
                as.annotation = TRUE)
}

# Numeric values (annotation=FALSE)
for(keyName in c('nCount_RNA', 'percent.mito', 'AvgDistanceToRef'))
{
    message('Adding ', keyName, '...')
    add_col_attr(loom, key = keyName,
                value = as.vector(seuratObj1[[keyName]][,1]),
                as.annotation = FALSE)
}

finalize(loom)
