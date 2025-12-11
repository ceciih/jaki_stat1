#!/bin/bash
#SBATCH --account=
#SBATCH --cluster=
#SBATCH --partition=
#SBATCH -t 1:00:00
#SBATCH --mem=500000
#SBATCH -o projection.out
# Load the Conda environment

module load cluster/wice/bigmem
source ./conda.sh
conda activate r40seurat40

R 

library(Seurat)
library(ggplot2)
library(patchwork)
reference <- readRDS("/$PATH/pbmc_multimodal_2023.rds")
p1 <- DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave(p1, file = "/lustre1/project/stg_00135/Cecilia/Data/umap.ref.png")
object <- readRDS('/lustre1/project/stg_00135/Cecilia/Data/p1ttmentannot.rds')
object <- SCTransform(object, verbose = TRUE)


anchors <- FindTransferAnchors(
  reference = reference,
  query = object,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

object <- MapQuery(
  anchorset = anchors,
  query = object,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

object <- TransferData(
  anchorset = anchors, 
  reference = reference,
  query = object,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT")
)
object <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference,
  query = object, 
  new.reduction.name = "ref.spca"
)
object <- ProjectUMAP(
  query = object, 
  query.reduction = "ref.spca", 
  reference = reference, 
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
ggsave(p1, file= "/lustre1/project/stg_00135/Cecilia/CITEseq/R_runs/level1.png")
ggsave(p2, file= "/lustre1/project/stg_00135/Cecilia/CITEseq/R_runs/level2.png")

saveRDS(object, file= '/lustre1/project/stg_00135/Cecilia/Data/projected.rds')
