#library(devtools)
#install_github('immunogenomics/presto')
#BiocManager::install("fgsea")
library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(dplyr)


#Prepare data ####
#Import Seurat object from basic Seurat pipeline
seuratObj <- readRDS("/lustre1/project/stg_00135/Cecilia/Data/ProjectedSatija/projected.both.rds")

dimall <- DimPlot(seuratObj, reduction="ref.umap", group.by="predicted.celltype.l1") + theme_void()
ggsave(dimall, file="/lustre1/project/stg_00135/Cecilia/Data/General_Study/test.png", width=5, height=5)

mono <- subset(seuratObj, subset = predicted.celltype.l1 == c("CD8 T"))
mono <- subset(mono, subset=subsample==c("a. Pre-Treatment P116RS", "b. 1 Week P116RS", "c. 4 Months P116RS"))

# Integrate ADT and RNA for analysis
DefaultAssay(mono) <- 'RNA'
mono <- NormalizeData(mono) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(mono) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(mono) <- rownames(mono[["ADT"]])
mono <- NormalizeData(mono, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

elvow <- ElbowPlot(mono, reduction = "pca") / ElbowPlot(mono, reduction = "apca")
ggsave(elvow, file='/lustre1/project/stg_00135/Cecilia/CITEseq/R_runs/elvow.png', width=21, height=7)

mono <- FindMultiModalNeighbors(
  mono, reduction.list = list("pca", "apca"), 
  dims.list = list(1:10, 1:7), modality.weight.name = "RNA.weight"
)

mono <- RunUMAP(mono, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mono <- FindClusters(mono, graph.name = "wsnn", algorithm = 3, resolution = 0.1, verbose = FALSE)

DefaultAssay(mono) <-"RNA"
#Creating artificial replicates in this data to use in analysis
mono$replicate <- c(rep("rep1",5257),rep("rep2",5257),
                       rep("rep3",5257),rep("rep4",5954),
                       rep("rep5",5954))
table(mono$subsample)

#Calculate DE stats of all genes ####
#Seurat Find(All)Markers doesn't provide stats for ALL genes
#Statistics via wilcoxauc() from presto package for KO versus WT
DEStats <- wilcoxauc(mono,group_by='subsample')
head(DEStats,2)

#Fill in the name of your favorite gene to check its results
DEStats[DEStats$feature=="STAT1",]

mdf <- msigdbr(species="Homo sapiens", category="C2",subcategory="CP:PID")
head(mdf)

#Split them into a list of gene sets
fgsea_sets <- split(mdf$gene_symbol,f=mdf$gs_name)

# Example for "1 Week P116RS"
DES_pre <- DEStats[DEStats$group == "a. Pre-Treatment P116RS", ]
DES_pre_sorted <- DES_pre[order(DES_pre$auc, decreasing = TRUE),]
ranks_pre <- DES_pre_sorted[,"auc"]
names(ranks_pre) <- DES_pre_sorted[,"feature"]

fgseaRes_pre <- fgsea(fgsea_sets, stats = ranks_pre)
# Example for "1 Week P116RS"
DES_1week <- DEStats[DEStats$group == "b. 1 Week P116RS", ]
DES_1week_sorted <- DES_1week[order(DES_1week$auc, decreasing = TRUE),]
ranks_1week <- DES_1week_sorted[,"auc"]
names(ranks_1week) <- DES_1week_sorted[,"feature"]

fgseaRes_1week <- fgsea(fgsea_sets, stats = ranks_1week)

# Repeat similar blocks for "4 Months P116RS"

# Example for "4 Week P116RS"
DES_4months <- DEStats[DEStats$group == "c. 4 Months P116RS", ]
DES_4months_sorted <- DES_4months[order(DES_4months$auc, decreasing = TRUE),]
ranks_4months <- DES_4months_sorted[,"auc"]
names(ranks_4months) <- DES_4months_sorted[,"feature"]

fgseaRes_4months <- fgsea(fgsea_sets, stats = ranks_4months)

# Example of how to compare the NES values across different timepoints
# Assuming fgseaRes_1week, fgseaRes_4months etc. are your GSEA results

compare_pathways_multiple <- function(pathway_names) {
  # Helper function to extract necessary statistics without leading-edge genes
  get_stats <- function(fgsea_res, pathway_name) {
    if (sum(fgsea_res$pathway == pathway_name) > 0) {
      entry <- fgsea_res[fgsea_res$pathway == pathway_name,]
      return(c(entry$NES, entry$pval, entry$padj, entry$size))
    } else {
      return(c(NA, NA, NA, NA))
    }
  }

  results <- list()
  
  # Process each pathway
  for (pathway in pathway_names) {
    pre_stats <- get_stats(fgseaRes_pre, pathway)
    week1_stats <- get_stats(fgseaRes_1week, pathway)
    months4_stats <- get_stats(fgseaRes_4months, pathway)
    
    # Create a data frame for each pathway
    results[[pathway]] <- data.frame(
      Timepoint = c("a. Pre-Treatment P116RS", "b. 1 Week P116RS", "c. 4 Months P116RS"),
      Pathway = pathway,
      NES = c(pre_stats[1], week1_stats[1], months4_stats[1]),
      P_Value = c(pre_stats[2], week1_stats[2], months4_stats[2]),
      Adjusted_P_Value = c(pre_stats[3], week1_stats[3], months4_stats[3]),
      Gene_Set_Size = c(pre_stats[4], week1_stats[4], months4_stats[4])
    )
  }
  # Combine all data frames into one
  combined_results <- do.call(rbind, results)
  return(combined_results)
}

#pathways_to_compare <- c("GOBP_T_HELPER_CELL_LINEAGE_COMMITMENT",
#"GOBP_POSITIVE_REGULATION_OF_T_HELPER_17_CELL_DIFFERENTIATION", 
#"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION",
#"GOBP_NEGATIVE_REGULATION_OF_T_HELPER_17_TYPE_IMMUNE_RESPONSE", 
#"GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION", 
#"GOBP_INTERLEUKIN_17_PRODUCTION",
#"GOBP_T_HELPER_17_CELL_DIFFERENTIATION",
#"GOBP_T_HELPER_17_CELL_LINEAGE_COMMITMENT",
#"GOBP_T_HELPER_17_TYPE_IMMUNE_RESPONSE",
#"GOBP_T_FOLLICULAR_HELPER_CELL_DIFFERENTIATION")

pathways_to_compare <- c(
  "GOBP_MONOCYTE_ACTIVATION",
  "GOBP_MONOCYTE_AGGREGATION",
  "GOBP_MONOCYTE_CHEMOTACTIC_PROTEIN_1_PRODUCTION",
  "GOBP_MONOCYTE_CHEMOTAXIS",
  "GOBP_MONOCYTE_DIFFERENTIATION",
  "GOBP_MONOCYTE_EXTRAVASATION")

pathways_to_compare <- c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_PI3k_AKT_MTOR_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

pathways_to_compare <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

# Monocyte relevant PID pathways (STAT1 GoF context)
pathways_to_compare <- c(
  "PID_IFNG_PATHWAY",
  "PID_IL1_PATHWAY",
  "PID_TNF_PATHWAY",
  "PID_NFKAPPAB_CANONICAL_PATHWAY",
  "PID_CASPASE_PATHWAY",
  "PID_FCER1_PATHWAY",
  "PID_LYSOPHOSPHOLIPID_PATHWAY",
  "PID_TOLL_ENDOGENOUS_PATHWAY",
  "PID_MYC_REPRESS_PATHWAY",
  "PID_IL6_7_PATHWAY",
  "PID_HIF1A_PATHWAY",
  "PID_P38_ALPHA_BETA_PATHWAY",
  "PID_MAPK_TRK_PATHWAY",
  "PID_P53_REGULATION_PATHWAY",
  "PID_NOTCH_PATHWAY",
  "PID_AP1_PATHWAY",
  "PID_FOXO_PATHWAY"
)

# T cell relevant PID pathways (STAT1 GoF context)
pathways_to_compare <- c(
  "PID_IFNG_PATHWAY",
  "PID_IL27_PATHWAY",
  "PID_IL2_STAT5_PATHWAY",
  "PID_TCR_PATHWAY",
  "PID_IL12_STAT4_PATHWAY",
  "PID_TCR_CALCIUM_PATHWAY",
  "PID_TCR_JNK_PATHWAY",
  "PID_TCR_RAS_PATHWAY",
  "PID_IL23_PATHWAY",
  "PID_IL4_2PATHWAY",
  "PID_IL5_PATHWAY",
  "PID_IL6_7_PATHWAY",
  "PID_FOXO_PATHWAY",
  "PID_P53_REGULATION_PATHWAY",
  "PID_NOTCH_PATHWAY",
  "PID_NFKAPPAB_CANONICAL_PATHWAY",
  "PID_AP1_PATHWAY",
  "PID_SMAD2_3PATHWAY",
  "PID_TGFBR_PATHWAY",
  "PID_ERBB_NETWORK_PATHWAY",
  "PID_PI3KCI_PATHWAY"
)
all_pathway_data <- compare_pathways_multiple(pathways_to_compare)
head(all_pathway_data)

library(ggplot2)

# Plot creation
gg <- ggplot(all_pathway_data, aes(x = Timepoint, y = Pathway, size = Gene_Set_Size, color = NES)) +
  geom_point(alpha=0.7) +  # Adjust transparency with alpha
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  scale_size(range = c(3, 12), name = "Gene Set Size") +
  labs(title = "Dot Plot of NES Across Timepoints and Pathways",
       x = "Timepoint",
       y = "Pathway") +
  theme_minimal(base_size = 14) +  # Ensures a minimal theme with a base font size of 14
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # Sets the panel background to white
    plot.background = element_rect(fill = "white", colour = "white"),  # Sets the plot background to white
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    axis.line = element_line(colour = "black"),  # Adds axis lines
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),  # Set x-axis labels at a 30-degree angle
    axis.ticks = element_line(color = "black")  # Ensures axis ticks are visible
  )

# Save the plot
ggsave(gg, file= "/lustre1/project/stg_00135/Cecilia/Data/General_Study/test.Mono.GO.png", width = 10, height = 10, units = "in")

pdf("/lustre1/project/stg_00135/Cecilia/Data/General_Study/Mono.GO.p1.pdf", width=8, height=6)
gg
dev.off()

########## P123RS_T385M

#library(devtools)
#install_github('immunogenomics/presto')
#BiocManager::install("fgsea")
library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(ggplot2)


#Prepare data ####
#Import Seurat object from basic Seurat pipeline
seuratObj <- readRDS("/lustre1/project/stg_00135/Cecilia/Data/projected.both.rds")

mono <- subset(seuratObj, subset = predicted.celltype.l1 == c("CD8 T"))
mono <- subset(mono, subset=subsample==c("d. Pre-Treatment P123RS", "e. 1 Week P123RS", "f. 4 Months P123RS"))

# Integrate ADT and RNA for analysis
DefaultAssay(mono) <- 'RNA'
mono <- NormalizeData(mono) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(mono) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(mono) <- rownames(mono[["ADT"]])
mono <- NormalizeData(mono, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


mono <- FindMultiModalNeighbors(
  mono, reduction.list = list("pca", "apca"), 
  dims.list = list(1:10, 1:7), modality.weight.name = "RNA.weight"
)

mono <- RunUMAP(mono, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mono <- FindClusters(mono, graph.name = "wsnn", algorithm = 3, resolution = 0.1, verbose = FALSE)

DefaultAssay(mono) <-"RNA"
#Creating artificial replicates in this data to use in analysis
mono$replicate <- c(rep("rep1",5257),rep("rep2",5257),
                       rep("rep3",5257),rep("rep4",5954),
                       rep("rep5",5954))
table(mono$subsample)

#Calculate DE stats of all genes ####
#Seurat Find(All)Markers doesn't provide stats for ALL genes
#Statistics via wilcoxauc() from presto package for KO versus WT
DEStats <- wilcoxauc(mono,group_by='subsample')
head(DEStats,2)

#Fill in the name of your favorite gene to check its results
DEStats[DEStats$feature=="STAT1",]

mdf <- msigdbr(species="Homo sapiens", category= "C2", subcategory="CP:PID")
head(mdf)

#Split them into a list of gene sets
fgsea_sets <- split(mdf$gene_symbol,f=mdf$gs_name)

# Example for "1 Week P116RS"
DES_pre <- DEStats[DEStats$group == "d. Pre-Treatment P123RS", ]
DES_pre_sorted <- DES_pre[order(DES_pre$auc, decreasing = TRUE),]
ranks_pre <- DES_pre_sorted[,"auc"]
names(ranks_pre) <- DES_pre_sorted[,"feature"]

fgseaRes_pre <- fgsea(fgsea_sets, stats = ranks_pre)
# Example for "1 Week P116RS"
DES_1week <- DEStats[DEStats$group == "e. 1 Week P123RS", ]
DES_1week_sorted <- DES_1week[order(DES_1week$auc, decreasing = TRUE),]
ranks_1week <- DES_1week_sorted[,"auc"]
names(ranks_1week) <- DES_1week_sorted[,"feature"]

fgseaRes_1week <- fgsea(fgsea_sets, stats = ranks_1week)

# Repeat similar blocks for "4 Months P116RS"

# Example for "4 Week P116RS"
DES_4months <- DEStats[DEStats$group == "f. 4 Months P123RS", ]
DES_4months_sorted <- DES_4months[order(DES_4months$auc, decreasing = TRUE),]
ranks_4months <- DES_4months_sorted[,"auc"]
names(ranks_4months) <- DES_4months_sorted[,"feature"]

fgseaRes_4months <- fgsea(fgsea_sets, stats = ranks_4months)

# Example of how to compare the NES values across different timepoints
# Assuming fgseaRes_1week, fgseaRes_4months etc. are your GSEA results

compare_pathways_multiple <- function(pathway_names) {
  # Helper function to extract necessary statistics without leading-edge genes
  get_stats <- function(fgsea_res, pathway_name) {
    if (sum(fgsea_res$pathway == pathway_name) > 0) {
      entry <- fgsea_res[fgsea_res$pathway == pathway_name,]
      return(c(entry$NES, entry$pval, entry$padj, entry$size))
    } else {
      return(c(NA, NA, NA, NA))
    }
  }

  results <- list()
  
  # Process each pathway
  for (pathway in pathway_names) {
    pre_stats <- get_stats(fgseaRes_pre, pathway)
    week1_stats <- get_stats(fgseaRes_1week, pathway)
    months4_stats <- get_stats(fgseaRes_4months, pathway)
    
    # Create a data frame for each pathway
    results[[pathway]] <- data.frame(
      Timepoint = c("d. Pre-Treatment P123RS", "e. 1 Week P123RS", "f. 4 Months P123RS"),
      Pathway = pathway,
      NES = c(pre_stats[1], week1_stats[1], months4_stats[1]),
      P_Value = c(pre_stats[2], week1_stats[2], months4_stats[2]),
      Adjusted_P_Value = c(pre_stats[3], week1_stats[3], months4_stats[3]),
      Gene_Set_Size = c(pre_stats[4], week1_stats[4], months4_stats[4])
    )
  }
  # Combine all data frames into one
  combined_results <- do.call(rbind, results)
  return(combined_results)
}

pathways_to_compare <- c("GOBP_T_HELPER_CELL_LINEAGE_COMMITMENT",
"GOBP_POSITIVE_REGULATION_OF_T_HELPER_17_CELL_DIFFERENTIATION", 
"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION",
"GOBP_NEGATIVE_REGULATION_OF_T_HELPER_17_TYPE_IMMUNE_RESPONSE", 
"GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION", 
"GOBP_INTERLEUKIN_17_PRODUCTION",
"GOBP_T_HELPER_17_CELL_DIFFERENTIATION",
"GOBP_T_HELPER_17_CELL_LINEAGE_COMMITMENT",
"GOBP_T_HELPER_17_TYPE_IMMUNE_RESPONSE",
"GOBP_T_FOLLICULAR_HELPER_CELL_DIFFERENTIATION")

pathways_to_compare <- c(
  "GOBP_MONOCYTE_ACTIVATION",
  "GOBP_MONOCYTE_AGGREGATION",
  "GOBP_MONOCYTE_CHEMOTACTIC_PROTEIN_1_PRODUCTION",
  "GOBP_MONOCYTE_CHEMOTAXIS",
  "GOBP_MONOCYTE_DIFFERENTIATION",
  "GOBP_MONOCYTE_EXTRAVASATION")

pathways_to_compare <- c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_PI3k_AKT_MTOR_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)


all_pathway_data <- compare_pathways_multiple(pathways_to_compare)


library(ggplot2)

# Plot creation
gg <- ggplot(all_pathway_data, aes(x = Timepoint, y = Pathway, size = Gene_Set_Size, color = NES)) +
  geom_point(alpha=0.7) +  # Adjust transparency with alpha
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  scale_size(range = c(3, 12), name = "Gene Set Size") +
  labs(title = "Dot Plot of NES Across Timepoints and Pathways",
       x = "Timepoint",
       y = "Pathway") +
  theme_minimal(base_size = 14) +  # Ensures a minimal theme with a base font size of 14
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # Sets the panel background to white
    plot.background = element_rect(fill = "white", colour = "white"),  # Sets the plot background to white
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    axis.line = element_line(colour = "black"),  # Adds axis lines
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),  # Set x-axis labels at a 30-degree angle
    axis.ticks = element_line(color = "black")  # Ensures axis ticks are visible
  )

# Save the plot
ggsave(gg, file= "/lustre1/project/stg_00135/Cecilia/Data/General_Study/test.Mono.GO.p2.png", width = 10, height = 10, units = "in")

pdf("/lustre1/project/stg_00135/Cecilia/Data/General_Study/Mono.GO.p2.pdf", width=8, height=6)
gg
dev.off()

### Combined Score

library(tidyverse)


# Calculate the combined score
data <- all_pathway_data %>%
  mutate(Combined_Score = NES * (1 / P_Value) * Gene_Set_Size)

# Calculate the average combined score for each base pathway
combined_scores <- data %>%
  group_by(rownames(data)) %>%
  summarise(Average_Combined_Score = mean(Combined_Score)) %>%
  arrange(desc(Average_Combined_Score))

# Create the bar plot
ggplot(combined_scores, aes(x = reorder(rownames(combined_scores), -Average_Combined_Score), y = Average_Combined_Score)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Average Combined Score for Each Pathway", x = "Pathway", y = "Average Combined Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()

# Show the plot
ggsave("/lustre1/project/stg_00135/Cecilia/Data/General/Study/average_combined_score_plot.mono.p1.pdf", width = 10, height = 6)