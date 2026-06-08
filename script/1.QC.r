library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
Con_1 <- readRDS("./Con_1.rds")
Con_2 <- readRDS("./Con_2.rds")
Con_3 <- readRDS("./Con_3.rds")
Con_4 <- readRDS("./Con_4.rds")
Con_5 <- readRDS("./Con_5.rds")
Con_6 <- readRDS("./Con_6.rds")
Cis_1 <- readRDS("./Cis_1.rds")
Cis_2 <- readRDS("./Cis_2.rds")
Cis_3 <- readRDS("./Cis_3.rds")
Cis_4 <- readRDS("./Cis_4.rds")
Cis_5 <- readRDS("./Cis_5.rds")
Cis_6 <- readRDS("./Cis_6.rds")
LRY15_1 <- readRDS("./LRY15_1.rds")
LRY15_2 <- readRDS("./LRY15_2.rds")
LRY15_3 <- readRDS("./LRY15_3.rds")
LRY15_4 <- readRDS("./LRY15_4.rds")
LRY15_5 <- readRDS("./LRY15_5.rds")
LRY15_6 <- readRDS("./LRY15_6.rds")
object_list <- list(
  Con_1   = Con_1,
  Con_2   = Con_2,
  Con_3   = Con_3,
  Con_4   = Con_4,
  Con_5   = Con_5,
  Con_6   = Con_6,
  Cis_1   = Cis_1,
  Cis_2   = Cis_2,
  Cis_3   = Cis_3,
  Cis_4   = Cis_4,
  Cis_5   = Cis_5,
  Cis_6   = Cis_6,
  LRY15_1 = LRY15_1,
  LRY15_2 = LRY15_2,
  LRY15_3 = LRY15_3,
  LRY15_4 = LRY15_4,
  LRY15_5 = LRY15_5,
  LRY15_6 = LRY15_6
)
sample_ids <- names(object_list)
for (i in seq_along(object_list)) {
  object_list[[i]]$sample <- sample_ids[i]
  object_list[[i]]$Batch <- sub("_\\d+$", "", sample_ids[i])
}
Bin50 <- merge(
  x = object_list[[1]],
  y = object_list[-1],
  add.cell.ids = sample_ids,
  project = "stereoseq"
)
Bin50$sample <- factor(
  Bin50$sample,
  levels = sample_ids
)
Bin50$Batch <- factor(
  Bin50$Batch,
  levels = c("Con", "Cis", "LRY15")
)
table(Bin50$sample)
table(Bin50$Batch)

#SCTransform normalization
DefaultAssay(Bin50) <- "RNA"
Bin50 <- SCTransform(
  object = Bin50,
  assay = "RNA",
  variable.features.n = 3000,
  verbose = FALSE
)
DefaultAssay(Bin50) <- "SCT"
Bin50 <- RunPCA(
  object = Bin50,
  assay = "SCT",
  npcs = 100,
  features = VariableFeatures(Bin50),
  verbose = FALSE
)
Bin50 <- RunHarmony(
  object = Bin50,
  group.by.vars = "sample",
  reduction = "pca",
  reduction.save = "harmony"
)
Bin50 <- FindNeighbors(
  object = Bin50,
  reduction = "harmony",
  dims = 1:20
)
Bin50 <- FindClusters(
  object = Bin50,
  resolution = 0.5
)
Bin50 <- RunUMAP(
  object = Bin50,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap"
)
spatial_matbatch <- Bin50@meta.data[,c("x","y")]
colnames(spatial_matbatch) <- c("spatial_1", "spatial_2")
spatial_objbatch <- CreateDimReducObject(
  embeddings = as.matrix(spatial_matbatch),
  key = "spatial_",
  assay = "SCT"
)
Bin50[["spatial"]] <- spatial_objbatch
saveRDS(Bin50,"./Bin50.rds")
#spatial gene number
dir.create("QC_spatial_plots", showWarnings = FALSE)
for (sample_name in names(object_list)) {
  
  p <- FeaturePlot(
    object_list[[sample_name]],
    features = "nFeature_RNA",
    reduction = "spatial"
  ) +
    scale_colour_gradientn(
      colours = c("#412169", "#fff200", "#ff0000"),
      breaks = c(1000, 1500, 2000, 2500)
    ) +
    ggtitle(sample_name)
  
  ggsave(
    filename = paste0("QC_spatial_plots/", sample_name, "_nFeature_RNA_spatial.pdf"),
    plot = p,
    width = 9,
    height = 8,
    units = "in"
  )
}
#gene_number_per_spot
gene_number <- FetchData(Bin50, vars = c("nFeature_RNA", "sample", "Batch"))
p_gene_number <- ggplot(
  gene_number,
  aes(x = sample, y = nFeature_RNA)
) +
  geom_boxplot(
    aes(color = factor(Batch)),
    fill = "white",
    outlier.shape = 16
  ) +
  labs(
    x = "Groups",
    y = "Gene number per spot"
  ) +
  scale_x_discrete(
    breaks = levels(gene_number$Batch)
  ) +
  scale_y_continuous() +
  scale_color_manual(
    values = c("#6CB8D2", "#479D88", "#E69F84")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    axis.title.x = element_text(
      hjust = 1,
      color = "black"
    ),
    axis.title.y = element_text(
      color = "black"
    ),
    axis.line = element_line(
      color = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    )
  )
ggsave(
  filename = "gene_number_per_spot.pdf",
  plot = p_gene_number,
  width = 5,
  height = 4,
  units = "in"
)
#mean_read_per_spot
read_perspot <- FetchData(Bin50, vars = c("nCount_RNA", "sample", "Batch"))
average_counts <- read_perspot %>% group_by(sample) %>% summarise(Average_ncount_RNA = mean(nCount_RNA))
average_counts$Batch <- c(rep("Con",6), rep("Cis",6), rep("LRY15", 6))
p_mean_read <- ggplot(
  average_counts,
  aes(
    x = sample,
    y = Average_ncount_RNA,
    fill = Batch
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  labs(
    x = "Group",
    y = "Mean read per spot"
  ) +
  scale_x_discrete(
    breaks = levels(average_counts$Batch)
  ) +
  scale_y_continuous() +
  scale_fill_manual(
    values = c("#479D88", "#6CB8D2", "#E69F84")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    axis.text.x = element_text(
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.line = element_line(
      color = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    )
  )
  ggsave(
  filename = "p_mean_read.pdf",
  plot = p_gene_read,
  width = 5,
  height = 4,
  units = "in"
)
#spot number
spot_number <- Bin50@meta.data %>%
  count(sample, Batch, name = "count")
spot_number$sample <- factor(
  spot_number$sample,
  levels = c(
    "Con_1", "Con_2", "Con_3", "Con_4", "Con_5", "Con_6",
    "Cis_1", "Cis_2", "Cis_3", "Cis_4", "Cis_5", "Cis_6",
    "LRY15_1", "LRY15_2", "LRY15_3", "LRY15_4", "LRY15_5", "LRY15_6"
  )
)
spot_number$Batch <- factor(
  spot_number$Batch,
  levels = c("Con", "Cis", "LRY15")
)
p_spot_number <- ggplot(
  spot_number,
  aes(
    x = sample,
    y = count,
    fill = Batch
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  labs(
    x = "Group",
    y = "Spot number"
  ) +
  scale_x_discrete(
    limits = levels(spot_number$sample)
  ) +
  scale_y_continuous() +
  scale_fill_manual(
    values = c(
      "Con" = "#479D88",
      "Cis" = "#6CB8D2",
      "LRY15" = "#E69F84"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      color = "black"
    ),
    axis.title.y = element_text(
      color = "black"
    ),
    axis.line = element_line(
      color = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    )
  )
ggsave(
  filename = "spot_number_per_sample.pdf",
  plot = p_spot_number,
  width = 8,
  height = 4,
  units = "in"
)
figures