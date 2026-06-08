library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(clusterProfiler)
library(org.Mm.eg.db)
Bin50 <- readRDS("./Bin50.rds")
all_markers_harmony_unannotation <- FindAllMarkers(Bin50, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, verbose = T)
all_markers_harmony_unannotation %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50_unannotation_harmony
write.table(top50_unannotation_harmony, file="./top50_unannotation_harmony.txt",sep = "\t")
new.cluster.ids <- c("Luteal cell", "Smooth muscle cell", "Theca-interstitial cell", "Epithelial cell", "Luteal cell", "Luteal cell", "Granulosa cell", "Stroma cell" ,"Granulosa cell", "Luteal cell", "Luteal cell", "Granulosa cell", "Oocyte", "Epithelial cell", "Oocyte")
names(new.cluster.ids) <- levels(Bin50)
names(new.cluster.ids) <- levels(Bin50@meta.data$seurat_clusters)
Bin50@meta.data$cell_type <- factor(Bin50@meta.data$seurat_clusters, 
                                             levels = names(new.cluster.ids), 
                                             labels = new.cluster.ids)
Idents(Bin50) <- Bin50$cell_type
#Gene set score for spot type
Oocyte_features <- list(c("Gdf9","Zp3","Ddx4", "Sycp3", "Ooep"))
Bin50 <- AddModuleScore(object = Bin50, features = Oocyte_features, ctrl = 100, name = "Oocyte_score")
data <- FetchData(Bin50, vars = c("Oocyte_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P1 <- ggplot(data, aes(x = Oocyte_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Oocyte_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Oocyte_score1, probs = c(0.01)), quantile(data$Oocyte_score1, probs = c(0.99))) + labs(x = "Oocyte marker gene set score", y = "Density")
ggsave(
  filename = "./P1.pdf",
  plot = P1,
  width = 5,
  height = 4,
  units = "in"
)
Granulosa_features <- list(c("Cyp19a1","Inhbb","Amh"))
Bin50 <- AddModuleScore(object = Bin50, features = Granulosa_features, ctrl = 100, name = "Granulosa_score")
data <- FetchData(Bin50, vars = c("Granulosa_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P2 <- ggplot(data, aes(x = Granulosa_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Granulosa_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Granulosa_score1, probs = c(0.01)), quantile(data$Granulosa_score1, probs = c(0.99))) + labs(x = "Granulosa marker gene set score", y = "Density")
ggsave(
  filename = "./P2.pdf",
  plot = P2,
  width = 5,
  height = 4,
  units = "in"
)
Theca_features <- list(c("Cyp17a1","Mgarp","Adh1"))
Bin50 <- AddModuleScore(object = Bin50, features = Theca_features, ctrl = 100, name = "Theca_score")
data <- FetchData(Bin50, vars = c("Theca_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P3 <- ggplot(data, aes(x = Theca_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Theca_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Theca_score1, probs = c(0.01)), quantile(data$Theca_score1, probs = c(0.99))) + labs(x = "Theca marker gene set score", y = "Density")
ggsave(
  filename = "./P3.pdf",
  plot = P3,
  width = 5,
  height = 4,
  units = "in"
)
Stroma_features <- list(c("Col1a1","Tcf21","Ptch1"))
Bin50 <- AddModuleScore(object = Bin50, features = Stroma_features, ctrl = 100, name = "Stroma_score")
data <- FetchData(Bin50, vars = c("Stroma_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P4 <- ggplot(data, aes(x = Stroma_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Stroma_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Stroma_score1, probs = c(0.01)), quantile(data$Stroma_score1, probs = c(0.99))) + labs(x = "Stroma marker gene set score", y = "Density")
ggsave(
  filename = "./P4.pdf",
  plot = P4,
  width = 5,
  height = 4,
  units = "in"
)
Smooth_muscle_cell_features <- list(c("Des","Acta2","Myh11"))
Bin50 <- AddModuleScore(object = Bin50, features = Smooth_muscle_cell_features, ctrl = 100, name = "Smooth_muscle_cell_score")
data <- FetchData(Bin50, vars = c("Smooth_muscle_cell_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P5 <- ggplot(data, aes(x = Smoothe_muscle_cell_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Smoothe_muscle_cell_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Smoothe_muscle_cell_score1, probs = c(0.01)), quantile(data$Smoothe_muscle_cell_score1, probs = c(0.99))) + labs(x = "Smoothe muscle cell marker gene set score", y = "Density")
ggsave(
  filename = "./P5.pdf",
  plot = P5,
  width = 5,
  height = 4,
  units = "in"
)
Luteal_cell_features <- list(c("Ptgfr","Neat1","Sfrp4"))
Bin50 <- AddModuleScore(object = Bin50, features = Luteal_cell_features, ctrl = 100, name = "Luteal_cell_score")
data <- FetchData(Bin50, vars = c("Luteal_cell_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P6 <- ggplot(data, aes(x = Luteal_cell_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Luteal_cell_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Luteal_cell_score1, probs = c(0.01)), quantile(data$Luteal_cell_score1, probs = c(0.99))) + labs(x = "Luteal cell marker gene set score", y = "Density")
ggsave(
  filename = "./P6.pdf",
  plot = P6,
  width = 5,
  height = 4,
  units = "in"
)
Immune_cell_features <- list(c("Cd74","Ptprc","Ctss"))
Bin50 <- AddModuleScore(object = Bin50, features = Immune_cell_features, ctrl = 100, name = "Immune_cell_score")
data <- FetchData(Bin50, vars = c("Immune_cell_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P7 <- ggplot(data, aes(x = Immune_cell_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Immune_cell_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Immune_cell_score1, probs = c(0.01)), quantile(data$Immune_cell_score1, probs = c(0.99))) + labs(x = "Immune cell marker gene set score", y = "Density")
ggsave(
  filename = "./P7.pdf",
  plot = P7,
  width = 5,
  height = 4,
  units = "in"
)
Epithelial_cell_features <- list(c("Krt19","Upk3b"))
Bin50 <- AddModuleScore(object = Bin50, features = Epithelial_cell_features, ctrl = 100, name = "Epithelial_cell_score")
data <- FetchData(Bin50, vars = c("Epithelial_cell_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P8 <- ggplot(data, aes(x = Epithelial_cell_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Epithelial_cell_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Epithelial_cell_score1, probs = c(0.01)), quantile(data$Epithelial_cell_score1, probs = c(0.99))) + labs(x = "Epithelial cell marker gene set score", y = "Density")
ggsave(
  filename = "./P8.pdf",
  plot = P8,
  width = 5,
  height = 4,
  units = "in"
)
Erythrocyte_features <- list(c("Hbb-bs","Hba-a1"))
Bin50 <- AddModuleScore(object = Bin50, features = Erythrocyte_features, ctrl = 100, name = "Erythrocyte_score")
data <- FetchData(Bin50, vars = c("Erythrocyte_score1", "cell_type"))
data$cell_type <- factor(data$cell_type, levels = c("Erythrocyte","Epithelial cell","Immune cell","Luteal cell","Smooth muscle cell","Stroma cell","Theca-interstitial cell","Granulosa cell","Oocyte"))
P9 <- ggplot(data, aes(x = Erythrocyte_score1, y = cell_type, fill = cell_type, color = cell_type)) + 
  geom_density_ridges(alpha = 0.7, size = 0.01) +
  geom_vline(xintercept = quantile(data$Erythrocyte_score1, probs = c(0.90)), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) + 
  scale_color_manual(values = c("Oocyte" = "#f8d90d", "Granulosa cell" = "#fc67fa", "Theca-interstitial cell" = "#00c6ff", "Stroma cell" = "#ff6581", "Smooth muscle cell" = "#b91d73", "Luteal cell" = "#c9d6ff", "Immune cell" = "#a5da66", "Epithelial cell" = "#ffd2d8", "Erythrocyte" = "#e578d6")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlim(quantile(data$Erythrocyte_score1, probs = c(0.01)), quantile(data$Erythrocyte_score1, probs = c(0.99))) + labs(x = "Erythrocyte marker gene set score", y = "Density")
ggsave(
  filename = "./P9.pdf",
  plot = P9,
  width = 5,
  height = 4,
  units = "in"
)
#GO enrichment for spot type
all_markers <- FindAllMarkers(Bin50, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, verbose = T)
all_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
clusters <- levels(top50$gene)
gene_lists <- lapply(clusters, function(cluster){top10_annotation$gene[top10_annotation$cluster ==cluster]})
GO_results <- lapply(gene_lists, function(genes){enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)})
write.table(GO_results, file = "./GO_results.txt", sep = "\t")
#spot type spatial distribution
colors <- c("#a8ff78","#f80759","#f7b733","#396afc","#ffc3a0","#f7ff00","#E0EAFC","#e74c3c","#ff5858","#182848","#ECE9E6","#5C258D","#5FC3E4","#16A085","#9796f0","#ec38bc")
p <- DimPlot(subset(Bin50, subset = group=='Con_1'), reduction = "spatial", cols = colors, pt.size = 0.2)
ggsave(
  filename = "./spot_type_spatial_distribution.pdf",
  plot = p,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Oocyte", "#f8d90d","#f7f7f8")
oocyte_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./oocyte_spatial.pdf",
  plot = oocyte_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Granulosa cell", "#fc67fa","#f7f7f8")
Granulosa_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Granulosa_cell_spatial.pdf",
  plot = Granulosa_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Theca-interstitial cell", "#00c6ff","#f7f7f8")
Theca_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Theca_cell_spatial.pdf",
  plot = Theca_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Stroma cell", "#ff6581","#f7f7f8")
Stroma_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Stroma_cell_spatial.pdf",
  plot = Stroma_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Smooth muscle cell", "#b91d73","#f7f7f8")
Smooth_muscle_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Smooth_muscle_cell_spatial.pdf",
  plot = Stroma_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Luteal cell", "#c9d6ff","#f7f7f8")
Luteal_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Luteal_cell_spatial.pdf",
  plot = Luteal_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Immune cell", "#a5da66","#f7f7f8")
Immune_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Immune_cell_spatial.pdf",
  plot = Immune_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Epithelial cell", "#ffd2d8","#f7f7f8")
Epithelial_cell_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Epithelial_cell_spatial.pdf",
  plot = Epithelial_cell_spatial,
  width = 5,
  height = 4,
  units = "in"
)
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$cell_type == "Erythrocyte", ""#e578d6","#f7f7f8")
Erythrocyte_spatial <- ggplot(
  Bin50@meta.data[Bin50@meta.data$sample == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = color_by_cluster
  )
) +
  geom_point(
    shape = 19,
    size = point_size
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed() +
  scale_color_identity()
ggsave(
  filename = "./Erythrocyte_spatial.pdf",
  plot = Erythrocyte_spatial,
  width = 5,
  height = 4,
  units = "in"
)
#Marker gene for spatial distribution
genes <- c("Gdf9", "Zp3")
groups <- c("Con", "Cis", "LRY15")

dir.create("FeaturePlot_PDF", showWarnings = FALSE)

for (gene in genes) {
  for (grp in groups) {
    
    p <- FeaturePlot(
      subset(Bin50, subset = Batch == grp),
      features = gene,
      raster = FALSE,
      reduction = "spatial"
    ) +
      scale_colour_gradientn(
        colours = c("#412169", "#a9d347", "#FDFC47")
      )
    
    ggsave(
      filename = paste0("FeaturePlot_PDF/", grp, "_", gene, ".pdf"),
      plot = p,
      width = 5,
      height = 4
    )
  }
}
genes <- c("Cyp19a1", "Inhbb")
groups <- c("Con", "Cis", "LRY15")

dir.create("FeaturePlot_PDF", showWarnings = FALSE)

for (gene in genes) {
  for (grp in groups) {
    
    p <- FeaturePlot(
      subset(Bin50, subset = Batch == grp),
      features = gene,
      raster = FALSE,
      reduction = "spatial"
    ) +
      scale_colour_gradientn(
        colours = c("#412169", "#a9d347", "#FDFC47")
      )
    
    ggsave(
      filename = paste0("FeaturePlot_PDF/", grp, "_", gene, ".pdf"),
      plot = p,
      width = 5,
      height = 4
    )
  }
}
genes <- c("Cyp17a1", "Mgarp")
groups <- c("Con", "Cis", "LRY15")

dir.create("FeaturePlot_PDF", showWarnings = FALSE)

for (gene in genes) {
  for (grp in groups) {
    
    p <- FeaturePlot(
      subset(Bin50, subset = Batch == grp),
      features = gene,
      raster = FALSE,
      reduction = "spatial"
    ) +
      scale_colour_gradientn(
        colours = c("#412169", "#a9d347", "#FDFC47")
      )
    
    ggsave(
      filename = paste0("FeaturePlot_PDF/", grp, "_", gene, ".pdf"),
      plot = p,
      width = 5,
      height = 4
    )
  }
}
genes <- c("Col1a1", "Tcf21")
groups <- c("Con", "Cis", "LRY15")

dir.create("FeaturePlot_PDF", showWarnings = FALSE)

for (gene in genes) {
  for (grp in groups) {
    
    p <- FeaturePlot(
      subset(Bin50, subset = Batch == grp),
      features = gene,
      raster = FALSE,
      reduction = "spatial"
    ) +
      scale_colour_gradientn(
        colours = c("#412169", "#a9d347", "#FDFC47")
      )
    
    ggsave(
      filename = paste0("FeaturePlot_PDF/", grp, "_", gene, ".pdf"),
      plot = p,
      width = 5,
      height = 4
    )
  }
}
genes <- c("Ptgfr", "Sfrp4")
groups <- c("Con", "Cis", "LRY15")

dir.create("FeaturePlot_PDF", showWarnings = FALSE)

for (gene in genes) {
  for (grp in groups) {
    
    p <- FeaturePlot(
      subset(Bin50, subset = Batch == grp),
      features = gene,
      raster = FALSE,
      reduction = "spatial"
    ) +
      scale_colour_gradientn(
        colours = c("#412169", "#a9d347", "#FDFC47")
      )
    
    ggsave(
      filename = paste0("FeaturePlot_PDF/", grp, "_", gene, ".pdf"),
      plot = p,
      width = 5,
      height = 4
    )
  }
}
saveRDS(Bin50,"./Bin50.rds")