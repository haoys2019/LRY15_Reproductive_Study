library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(grid)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(CytoTRACE)
library(dorothea)
library(viper)
library(pheatmap)
Bin50 <- readRDS("./Bin50.rds")
#GO enrichment
Bin50@meta.data$group_type <- paste(Bin50@meta.data$Batch, Bin50_merge@meta.data$cell_type, sep = "_")
head(Bin50_merge@meta.data)
Cis_Con_granulosa_cell_marker <- FindMarkers(Bin50, group.by = "group_type", ident.1 = "Cis_Granulosa cell", ident.2 = "Con_Granulosa cell", min.pct = 0.25, logfc.threshold = 0.25)
LRY15_Cis_granulosa_cell_marker <- FindMarkers(Bin50, group.by = "group_type", ident.1 = "LRY15_Granulosa cell", ident.2 = "Cis_Granulosa cell", min.pct = 0.25, logfc.threshold = 0.25)
Cis_Con_granulosa_cell_marker$gene <- rownames(granulosa_cell_marker)
CiS_Con_granulosa_cell_up <- Cis_Con_granulosa_cell_marker[
  Cis_Con_granulosa_cell_marker$avg_log2FC > 0.25 &
    Cis_Con_granulosa_cell_marker$p_val_adj < 0.05,
]
CiS_Con_granulosa_cell_down <- Cis_Con_granulosa_cell_marker[
  Cis_Con_granulosa_cell_marker$avg_log2FC < 0.25 &
    Cis_Con_granulosa_cell_marker$p_val_adj < 0.05,
]
CiS_CON_DEGs <- Cis_Con_granulosa_cell_marker[
  abs(Cis_Con_granulosa_cell_marker$avg_log2FC) > 0.25 &
    Cis_Con_granulosa_cell_marker$p_val_adj < 0.05,
]
LRY15_Cis_granulosa_cell_marker$gene <- rownames(LRY15_Cis_granulosa_cell_marker)
LRY15_Cis_granulosa_cell_up <- LRY15_Cis_granulosa_cell_marker[
  LRY15_Cis_granulosa_cell_marker$avg_log2FC > 0.25 &
    LRY15_Cis_granulosa_cell_marker$p_val_adj < 0.05,
]
LRY15_Cis_granulosa_cell_down <- LRY15_Cis_granulosa_cell_marker[
  LRY15_Cis_granulosa_cell_marker$avg_log2FC < 0.25 &
    LRY15_Cis_granulosa_cell_marker$p_val_adj < 0.05,
]
LRY15_Cis_DEGs <- LRY15_Cis_granulosa_cell_marker[
  abs(LRY15_Cis_granulosa_cell_marker$avg_log2FC) > 0.25 &
    LRY15_Cis_granulosa_cell_marker$p_val_adj < 0.05,
]
CiS_Con_granulosa_cell_down_v <- rownames(CiS_Con_granulosa_cell_down)
LRY15_Cis_granulosa_cell_up_v <- rownames(LRY15_Cis_granulosa_cell_up)
venn_plot_granulosa_cell_Cis_down_LRY15_up <- draw.pairwise.venn(
  area1 = length(CiS_Con_granulosa_cell_down_v),
  area2 = length(LRY15_Cis_granulosa_cell_up_v),
  cross.area = length(intersect(CiS_Con_granulosa_cell_down_v, LRY15_Cis_granulosa_cell_up_v)),
  category = c("CIS vs CON down", "LRY15 vs CIS up"),
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("lightblue", "#f0cb35")
)
pdf("Venn_CIS_down_LRY15_up.pdf", width = 5, height = 5)
grid.draw(venn_plot_granulosa_cell_Cis_down_LRY15_up)
dev.off()
CiS_Con_granulosa_cell_up_v <- rownames(CiS_Con_granulosa_cell_up)
LRY15_Cis_granulosa_cell_down_v <- rownames(LRY15_Cis_granulosa_cell_down)
venn_plot_granulosa_cell_Cis_up_LRY15_down <- draw.pairwise.venn(
  area1 = length(CiS_Con_granulosa_cell_up_v),
  area2 = length(LRY15_Cis_granulosa_cell_down_v),
  cross.area = length(intersect(CiS_Con_granulosa_cell_up_v, LRY15_Cis_granulosa_cell_down_v)),
  category = c("CIS vs CON up", "LRY15 vs CIS down"),
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("lightblue", "#f0cb35")
)
pdf("Venn_CIS_up_LRY15_dowm.pdf", width = 5, height = 5)
grid.draw(venn_plot_granulosa_cell_Cis_up_LRY15_down)
dev.off()
ego_CiS_Con_granulosa_cell_down <- clusterProfiler::enrichGO(gene = CiS_Con_granulosa_cell_down_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_CiS_Con_granulosa_cell_up <- clusterProfiler::enrichGO(gene = CiS_Con_granulosa_cell_up_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_LRY15_Cis_granulosa_cell_up <- clusterProfiler::enrichGO(gene = LRY15_Cis_granulosa_cell_up_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_LRY15_Cis_granulosa_cell_down <- clusterProfiler::enrichGO(gene = LRY15_Cis_granulosa_cell_down_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
x <- as.data.frame(ego_LRY15_Cis_granulosa_cell_down_specfic_pathway)
x$log10_pvalue <- -log10(x$pvalue)
x$Description <- factor(
  x$Description,
  levels = rev(unique(x$Description))
)
p_GO <- ggplot(
  x,
  aes(x = log10_pvalue, y = Description)
) +
  geom_segment(
    aes(
      x = 0,
      xend = log10_pvalue,
      y = Description,
      yend = Description
    ),
    color = "grey70",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      size = log10_pvalue,
      color = Count
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    name = "-Log10(P-value)",
    range = c(3, 10)
  ) +
  scale_color_gradient(
    low = "#a5b9e9",
    high = "#1307ff",
    name = "Gene Count"
  ) +
  labs(
    x = "-Log10(P-value)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  )
ggsave(
  filename = "GO__LRY15_Cis_granulosa_cell_down_specific_pathway.pdf",
  plot = p_GO,
  width = 7,
  height = 6,
  units = "in"
)
#Gene set score for ovulation
Granulosa_cell <- subset(Bin50, subset = cell_type == "Granulosa cell")
Ovulation <- read_excel("Ovulation.xlsx")
gene <- as.list(Ovulation)
Granulosa_cell <- AddModuleScore(object = Granulosa_cell, features = gene, ctrl = 100, name = 'ovulation')
colnames(Granulosa_cell@meta.data)
Idents(Granulosa_cell) <- Granulosa_cell$Batch
data <- FetchData(Granulosa_cell, vars = c("Batch", "Ovulation1"))
test_Con_Cis <- wilcox.test(Ovulation1 ~ Batch, data = data, subset = Batch %in% c("Con", "Cis"))
test_LRY15_Cis <- wilcox.test(Ovulation1 ~ Batch, data = data, subset = Batch %in% c("LRY15", "Cis"))
test_Con_LRY15 <- wilcox.test(Ovulation1 ~ Batch, data = data, subset = Batch %in% c("LRY15", "Con"))
ggplot(data, aes(Batch, Ovulation1, color = Batch))+geom_boxplot()+theme_bw()+RotatedAxis() + geom_text(aes(x = 1.5, y = 0.1, label = paste("p = ", format.pval(test_Con_Cis$p.value, digits = 3))), vjust = -0.5)+geom_text(aes(x = 2.5, y = 0.1, label = paste("p = ", format.pval(test_LRY15_Cis$p.value, digits = 3))), vjust = -0.5)+geom_text(aes(x = 2, y = 0.18, label = paste("p = ", format.pval(test_Con_LRY15$p.value, digits = 3))), vjust = -0.5) + scale_color_manual(values = c("#D55640","#6CB8D2","#E69F84"))+labs(x = "Granulosa cell", y = "Ovulation1")+theme(legend.background = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylim(-0.1, 0.3)
#Spatial gene set score 
Con_1_granulosa_cell <- subset(Grabulosa_cell, subset = sample == "Con_1")
data <- FetchData(Con_1_granulosa_cell, vars = c("x", "y", "Ovulation1"))
p_ovulation <- ggplot() +
  geom_point(
    data = Con_1_granulosa_cell@meta.data[Con_1_granulosa_cell@meta.data$sample == "Con_1", ],
    aes(
      x = x,
      y = y
    ),
    color = "#f7f6f5",
    size = 1
  ) +
  geom_point(
    data = data,
    aes(
      x = x,
      y = y,
      color = Ovulation1
    ),
    size = 0.2
  ) +
  scale_color_gradientn(
    colors = c(
      "#0052d4",
      "#4364f7",
      "#6fb1fc",
      "yellow",
      "orange",
      "red"
    ),
    breaks = c(-0.4,-0.2, 0, 0.2, 0.4),
    limits = c(-0.4, 0.4),
    name = "Score"
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right"
  )
ggsave(
  filename = "Con_1_ovulation_score_spatial.pdf",
  plot = p_ovulation,
  width = 5,
  height = 4,
  units = "in"
)
#spatial gene expression
data <- FetchData(Con1_granulosa_cell, vars = c("x", "y", "Inhba"))
p_Inhba_spatial <- ggplot() + geom_point(data = Con_1@meta.data[Con_1@meta.data$group == "Con_1", ], aes(x = x, y =y), color= "#f7f7f8", size =1) + geom_point(data = data, aes(x = x, y = y, color = Inhba), size = 0.2) + scale_color_gradientn(colors = c("#0575e6", "cyan", "green", "yellow", "orange","red" ),breaks =c(0, 1, 2, 3), limits = c(0,3)) + coord_fixed() + theme_void() + theme(legend.position = "right") + labs(color = "Inhba")
ggsave(
  filename = "Inhba_spatial.pdf",
  plot = p_Inhba_spatial,
  width = 5,
  height = 4,
  units = "in"
)
test_Con_Cis <- wilcox.test(Inhba ~ Batch, data = data, subset = Batch %in% c("Con", "Cis"))
test_LRY15_Cis <- wilcox.test(Inhba ~ Batch, data = data, subset = Batch %in% c("LRY15", "Cis"))
test_Con_LRY15 <- wilcox.test(Inhba ~ Batch, data = data, subset = Batch %in% c("LRY15", "Con"))
aggregate(Inhba ~ Batch, data = data, mean)
p_Inhba_violin <- ggplot(
  data,
  aes(
    x = Batch,
    y = Inhba,
    fill = Batch,
    color = Batch
  )
) +
  geom_violin() +
  scale_fill_manual(
    values = c("#D55640", "#6CB8D2", "#E69F84")
  ) +
  scale_color_manual(
    values = c("#D55640", "#6CB8D2", "#E69F84")
  ) +
  labs(
    x = "",
    y = "Expression level",
    title = "Inhba"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      color = "black"
    ),
    axis.text.y = element_text(
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
    ),
    legend.position = "none"
  )
ggsave(
  filename = "Inhba_violin.pdf",
  plot = p_Inhba_violin,
  width = 5,
  height = 4,
  units = "in"
)
#Multiple violin plot 
p_violin <- ggplot(data, aes(x = Batch, y = Gja1, fill = Batch)) +
  geom_violin(width = 0.5, trim = TRUE, adjust = 1.3, bw = 0.35, color = NA) +
  stat_summary(fun = median, geom = "point", size = 2, color = "white") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.05, color = "black", coef = Inf, linewidth = 0.3) +
  scale_fill_manual(values = c("#D55640", "#6CB8D2", "#E69F84")) +
  labs(x = NULL, y = "Expression level", title = "Gja1") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )
lims <- c(1,4)#y轴坐标
brks <- seq(-1:4, by = 1)
rmY <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
P1 <- P1 + scale_y_continuous(limits = lims, breaks = brks)
P2 <- P2 + scale_y_continuous(limits = lims, breaks = brks) + rmY
P3 <- P3 + scale_y_continuous(limits = lims, breaks = brks) + rmY
P4 <- P4 + scale_y_continuous(limits = lims, breaks = brks) + rmY
(P1 | P2 | P3 | P4) + plot_layout(ncol = 4, guides = "collect")
#Granulosa cell subtype
Granulosa_cell <- subset(Bin50, subset = cell_type == "Granulosa cell")
DefaultAssay(Granulosa_cell) <- "RNA"
Granulosa_cell <- SCTransform(object = Granulosa_cell, assay = "RNA", variable.features.n = 3000,verbose = FALSE)
DefaultAssay(Granulosa_cell) <- "SCT"
Granulosa_cell <- RunPCA(object = Granulosa_cell, assay = "SCT", npcs = 50, features = VariableFeatures(Granulosa_cell), verbose = FALSE)
Granulosa_cell <- RunHarmony(object = Granulosa_cell, group.by.vars = "sample", reduction = "pca", reduction.save = "harmony")
Granulosa_cell <- FindNeighbors(object = Granulosa_cell, reduction = "harmony", dims = 1:20)
Granulosa_cell <- FindClusters(object = Granulosa_cell, resolution = 0.2)
Granulosa_cell <- RunUMAP(object = Granulosa_cell, reduction = "harmony", dims = 1:20, reduction.name = "umap")
sub_Granulosa_cell_all_markers <- FindAllMarkers(Granulosa_cell, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5, verbose = T)
list_subgranulosa_cell <- split(sub_Granulosa_cell_all_markers$gene, sub_Granulosa_cell_all_markers$cluster)
granulosa_1_gene <- list_subgranulosa_cell$Granulosa_1
granulosa_2_gene <- list_subgranulosa_cell$Granulosa_2
granulosa_3_gene <- list_subgranulosa_cell$Granulosa_3
granulosa_4_gene <- list_subgranulosa_cell$Granulosa_4
ego_granulosa_1 <- clusterProfiler::enrichGO(
  gene          = granulosa_1_gene,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "SYMBOL",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ego_granulosa_2 <- clusterProfiler::enrichGO(
  gene          = granulosa_1_gene,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "SYMBOL",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ego_granulosa_3 <- clusterProfiler::enrichGO(
  gene          = granulosa_1_gene,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "SYMBOL",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ego_granulosa_4 <- clusterProfiler::enrichGO(
  gene          = granulosa_1_gene,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "SYMBOL",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
p_heatmap_sub_granulosa <- DoHeatmap(
  object = sub_Granulosa_cell,
  features = all_marker_sub_Granulosa_cell$gene,
  label = FALSE
) +
  scale_fill_gradient2(
    low = "#3f2b96",
    mid = "#FFFFFF",
    high = "red"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )
#Granulosa cell subtype spatial distribution
Con_1@meta.data$granulosacell <- ifelse(Con_1@meta.data$granulosacell == "Granulosa cell", "#f7f7f8", "#f7f7f8")
p <- ggplot(
  data = Con_1@meta.data[Con_1@meta.data$group == "Con_1", ],
  aes(
    x = x,
    y = y,
    color = granulosacell
  )
) +
  geom_point(
    shape = 19,
    size = 1
  ) +
  scale_color_identity() +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "none"
  )
sub_Granulosa_cell@meta.data$color_by_celltype <- ifelse(
  sub_Granulosa_cell@meta.data$Granulosa_type == "Granulosa_1",
  "#ffffd2",
  
  ifelse(
    sub_Granulosa_cell@meta.data$Granulosa_type == "Granulosa_2",
    "#fcbad3",
    
    ifelse(
      sub_Granulosa_cell@meta.data$Granulosa_type == "Granulosa_3",
      "#aa96da",
      
      ifelse(
        sub_Granulosa_cell@meta.data$Granulosa_type == "Granulosa_4",
        "#a8d8ea",
        "#F7F7F8"
      )
    )
  )
)
p_granulosa_subtype <- p +
  geom_point(
    data = sub_Granulosa_cell@meta.data[
      sub_Granulosa_cell@meta.data$group == "Con_1",
    ],
    aes(
      x = x,
      y = y,
      color = color_by_celltype
    ),
    shape = 19,
    size = 0.2
  ) +
  scale_color_identity()
ggsave(
  filename = "Con_1_granulosa_subtype_spatial.pdf",
  plot = p_granulosa_subtype,
  width = 5,
  height = 4,
  units = "in"
)
saveRDS(Granulosa_cell, file = "./Granulosa_cell.rds")
#CytoTRACE
exp1 <- as.matrix(Granulosa_cell@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- Granulosa_cell$clusters
phenot <- as.character(phenot)
names(phenot) <- rownames(Granulosa_cell@meta.data)
emb <- Granulosa_cell@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')
#DoRothEA
sub_Granulosa_cell <- readRDS("~/_Granulosa_cell.rds")
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
regulon <- dorothea_regulon_mouse %>% dplyr::filter(confidence %in% c("A", "B", "C"))
sub_Granulosa_cell <- run_viper(sub_Granulosa_cell, regulon, options = list(method = "scale", minsize = 4, eset.filter = FALSE, cores = 5, verbose = FALSE))
Assays(sub_Granulosa_cell)
sub_Granulosa_cell@assays$dorothea@data[1:4,1:4]
DefaultAssay(object = sub_Granulosa_cell) <- "dorothea"
table(Idents(sub_Granulosa_cell))
sub_Granulosa_cell <- ScaleData(sub_Granulosa_cell)
viper_scores_df <- GetAssayData(sub_Granulosa_cell, layer = "scale.data", assay = "dorothea") %>% data.frame(check.names = F) %>% t()
CellsClusters <- data.frame(cell = names(Idents(sub_Granulosa_cell)), cell_type = as.character(sub_Granulosa_cell@meta.data$group_granulosa_type), check.names = F)
View(CellsClusters)
sub_Granulosa_cell@meta.data$group_granulosa_type <- factor(sub_Granulosa_cell@meta.data$group_granulosa_type, levels = c("Con_Granulosa_1", "Cis_Granulosa_1", "LRB5_Granulosa_1", "Con_Granulosa_2", "Cis_Granulosa_2", "LRB5_Granulosa_2", "Con_Granulosa_3", "Cis_Granulosa_3", "LRB5_Granulosa_3", "Con_Granulosa_4", "Cis_Granulosa_4", "LRB5_Granulosa_4"))
CellsClusters <- data.frame(cell = names(Idents(sub_Granulosa_cell)), cell_type = as.character(sub_Granulosa_cell@meta.data$group_granulosa_type), check.names = F)
CellsClusters <- CellsClusters[order(CellsClusters$cell_type),]
View(CellsClusters)
viper_scores_clusters <- viper_scores_df %>% data.frame() %>% rownames_to_column("cell") %>% gather(tf, activity, -cell) %>% inner_join(CellsClusters)
summarized_viper_scores <- viper_scores_clusters %>% group_by(tf, cell_type) %>% summarise(avg = mean(activity), std = sd(activity))
View(summarized_viper_scores)
highly_variable_tfs <- summarized_viper_scores %>% group_by(tf) %>% mutate(var = var(avg)) %>% ungroup() %>% top_n(240, var) %>% distinct(tf)
View(highly_variable_tfs)
summarized_viper_scores_df <- summarized_viper_scores %>% semi_join(highly_variable_tfs, by = "tf") %>% dplyr::select(-std) %>% spread(tf,avg) %>% data.frame(row.names = 1, check.names = FALSE)
View(summarized_viper_scores_df)
summarized_viper_scores_df[1:4,1:4]
View(summarized_viper_scores)
palette_length = 100
my_color = colorRampPalette(c("navy", "white","red"))(palette_length)
my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
t_summarized_viper_scores_df <- t(summarized_viper_scores_df)
p_dorothea <- pheatmap(t_summarized_viper_scores_df,fontsize = 14, fontsize_row = 10, color = my_color, breaks = my_breaks, main = "DoRothEA(ABC)", angle_col = 45, treeheight_col = 0, broder_color =NA,cluster_cols = FALSE)
ggsave(
  filename = "heatmap_dorothea.pdf",
  plot = p_dorothea,
  width = 5,
  height = 4,
  units = "in"
)