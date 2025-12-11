library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

object <- readRDS("/lustre1/project/stg_00135/Cecilia/Data/projected.both.rds")

# Assuming `reference_dataset@reductions$pca@cell.embeddings` exists
reference <- readRDS("/lustre1/project/stg_00135/Cecilia/Data/pbmc_multimodal_2023.rds")
reference <- subset(x = reference, subset = donor == c("P3", "P6", "P7"))
reference@reductions
# Extract PCA embeddings
query_embeddings <- Embeddings(object, "ref.spca")
ref_embeddings <- Embeddings(reference, "spca")

# Function to calculate average distance to reference cells
calculate_average_distances <- function(query_embedding, ref_embeddings) {
  distances <- apply(ref_embeddings, 1, function(ref) sqrt(sum((query_embedding - ref)^2)))
  mean(distances)
}

# Apply this function to each cell in the query dataset
# Run this in cluster, takes long: average_distances <- apply(query_embeddings, 1, calculate_average_distances, ref_embeddings = ref_embeddings)

average_distances <- read.csv2("/lustre1/project/stg_00135/Cecilia/Data/average.dist.csv")
# Split the 'X.x' column on the comma
split_data <- strsplit(as.character(average_distances$X.x), ",")
# Create a new data frame with two columns: 'cell_id' and 'distance'
average_distances <- data.frame(
  cell_id = sapply(split_data, `[`, 1),
  distance = as.numeric(sapply(split_data, `[`, 2))
)

# Check the resulting data frame
head(average_distances)

# Add the distance data as metadata. Ensure the row names of the object match those in 'cell_id'.
rownames(average_distances) <- average_distances$cell_id
object <- AddMetaData(object, metadata = average_distances$distance, col.name = "AvgDistanceToRef")

# Check if UMAP is available, if not, compute it
if (!"ref.umap" %in% names(object@reductions)) {
  object <- RunUMAP(object, dims = 1:30)
}

# Plot using FeaturePlot
dist <- FeaturePlot(object, features = "AvgDistanceToRef", reduction="ref.umap", split.by="subsample") +scale_color_viridis()
ggsave(dist, file= "/lustre1/project/stg_00135/Cecilia/CITEseq/R_runs/dist2ref.both.png", width=30)



##### Adding from DISCO
table(object$subsample)
object1 <- subset(x = object, subset = subsample == c("a. Pre-Treatment P116RS"))
object2 <- subset(x = object, subset = subsample == c("d. Pre-Treatment P123RS"))

rm(object)

# P116RS
cellmapper_res1 = read.csv("/lustre1/project/stg_00135/Cecilia/Data/DISCO_Dist2Ref_P116RS-Pre.txt", sep = "\t")
rownames(cellmapper_res1) = cellmapper_res1$cell
object1 = AddMetaData(object1, cellmapper_res1) #rna is your input Seurat object
dim <- FeaturePlot(object1, features="cellmapper_distance_to_reference", max.cutoff=20)
dim <- dim + scale_color_viridis() + theme_void()
pdf("/lustre1/project/stg_00135/Cecilia/Data/General_Study/Dist2Ref_DISCO_P116RS_Pre.pdf",width=5, height=5)
dim
dev.off()
# P123RS
cellmapper_res2 = read.csv("/lustre1/project/stg_00135/Cecilia/Data/DISCO_Dist2Ref_P123RS-Pre.txt", sep = "\t")
rownames(cellmapper_res2) = cellmapper_res2$cell
object2 = AddMetaData(object2, cellmapper_res2) #rna is your input Seurat object
dim2 <- FeaturePlot(object2, features="cellmapper_distance_to_reference", max.cutoff=20)
dim2 <- dim2 + scale_color_viridis() + theme_void()

pdf("/lustre1/project/stg_00135/Cecilia/Data/General_Study/Dist2Ref_DISCO_P123RS_Pre.pdf",width=5, height=5)
dim2
dev.off()

dim3 <- DimPlot(object, group.by="predicted.celltype.l1", label=T) + theme_void() + NoLegend()
ggsave(dim3, file="/lustre1/project/stg_00135/Cecilia/CITEseq/R_runs/DimUmap_All.png", width=5, height=5)
colnames(object1@meta.data)
