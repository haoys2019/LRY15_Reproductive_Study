library(Seurat)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(CellChat)
Bin50 <- readRDS("./Bin50.rds")
#GO enrichment
Bin50@meta.data$group_type <- paste(Bin50@meta.data$Batch, Bin50_merge@meta.data$cell_type, sep = "_")
head(Bin50_merge@meta.data)
Cis_Con_oocytemarker <- FindMarkers(Bin50, group.by = "group_type", ident.1 = "Cis_Oocyte", ident.2 = "Con_Oocyte", min.pct = 0.25, logfc.threshold = 0.25)
LRY15_Cis_oocytemarker <- FindMarkers(Bin50, group.by = "group_type", ident.1 = "LRY15_Oocyte", ident.2 = "Cis_Oocyte", min.pct = 0.25, logfc.threshold = 0.25)
Cis_Con_oocytemarker$gene <- rownames(Cis_Con_oocytemarker)
CiS_Con_oocyte_up <- Cis_Con_oocytemarker[
  Cis_Con_oocytemarker$avg_log2FC > 0.25 &
    Cis_Con_oocytemarker$p_val_adj < 0.05,
]
CiS_Con_oocyte_down <- Cis_Con_oocytemarker[
  Cis_Con_oocytemarker$avg_log2FC < -0.25 &
    Cis_Con_oocytemarker$p_val_adj < 0.05,
]
CiS_CON_DEGs <- Cis_Con_oocytemarker[
  abs(Cis_Con_oocytemarker$avg_log2FC) > 0.25 &
    Cis_Con_oocytemarker$p_val_adj < 0.05,
]
LRY15_Cis_oocytemarker$gene <- rownames(LRY15_Cis_oocytemarker)
LRY15_Cis_oocyte_up <- LRY15_Cis_oocytemarker[
  LRY15_Cis_oocytemarker$avg_log2FC > 0.25 &
    LRY15_Cis_oocytemarker$p_val_adj < 0.05,
]
LRY15_Cis_oocyte_down <- LRY15_Cis_oocytemarker[
  LRY15_Cis_oocytemarker$avg_log2FC < -0.25 &
    LRY15_Cis_oocytemarker$p_val_adj < 0.05,
]
LRY15_Cis_DEGs <- LRY15_Cis_oocytemarker[
  abs(LRY15_Cis_oocytemarker$avg_log2FC) > 0.25 &
    LRY15_Cis_oocytemarker$p_val_adj < 0.05,
]
CiS_Con_oocyte_down_v <- rownames(CiS_Con_oocyte_down)
LRY15_Cis_oocyte_up_v <- rownames(LRY15_Cis_oocyte_up)
venn_plot_oocyte_Cis_down_LRY15_up <- draw.pairwise.venn(
  area1 = length(CiS_Con_oocyte_down_v),
  area2 = length(LRY15_Cis_oocyte_up_v),
  cross.area = length(intersect(CiS_Con_oocyte_down_v, LRY15_Cis_oocyte_up_v)),
  category = c("CIS vs CON down", "LRY15 vs CIS up"),
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("lightblue", "#f0cb35")
)
pdf("Venn_CIS_down_LRY15_up.pdf", width = 5, height = 5)
grid.draw(venn_plot_oocyte_Cis_down_LRY15_up)
dev.off()
CiS_Con_oocyte_up_v <- rownames(CiS_Con_oocyte_up)
LRY15_Cis_oocyte_down_v <- rownames(LRY15_Cis_oocyte_down)
venn_plot_oocyte_Cis_up_LRY15_down <- draw.pairwise.venn(
  area1 = length(CiS_Con_oocyte_up_v),
  area2 = length(LRY15_Cis_oocyte_down_v),
  cross.area = length(intersect(CiS_Con_oocyte_up_v, LRY15_Cis_oocyte_down_v)),
  category = c("CIS vs CON up", "LRY15 vs CIS down"),
  fill = c("red", "blue"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("lightblue", "#f0cb35")
)
pdf("Venn_CIS_up_LRY15_dowm.pdf", width = 5, height = 5)
grid.draw(venn_plot_oocyte_Cis_up_LRY15_down)
dev.off()
ego_CiS_Con_oocyte_down <- clusterProfiler::enrichGO(gene = CiS_Con_oocyte_down_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_CiS_Con_oocyte_up <- clusterProfiler::enrichGO(gene = CiS_Con_oocyte_up_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_LRY15_Cis_oocyte_up <- clusterProfiler::enrichGO(gene = LRY15_Cis_oocyte_up_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
ego_LRY15_Cis_oocyte_down <- clusterProfiler::enrichGO(gene = LRY15_Cis_oocyte_down_v,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
x <- as.data.frame(ego_LRY15_Cis_oocyte_down_specfic_pathway)
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
  filename = "GO_lollipop_LRY15_Cis_oocyte_down_specific_pathway.pdf",
  plot = p_GO_lollipop,
  width = 7,
  height = 6,
  units = "in"
)
#Gene set score
oocyte <- subset(Bin50, subset = cell_type == "Oocyte")
intrinsic_apoptosis_by_P53 <- read_excel("intrinsic_apoptosis_by_P53.xlsx")
gene <- as.list(intrinsic_apoptosis_by_P53)
oocyte <- AddModuleScore(object = oocyte, features = gene, ctrl = 100, name = 'intrinsic_apoptosis_by_P53')
colnames(oocyte@meta.data)
Idents(oocyte) <- oocyte$Batch
data <- FetchData(oocyte, vars = c("Batch", "intrinsic_apoptosis_by_P531"))
test_Con_Cis <- wilcox.test(intrinsic_apoptosis_by_P531 ~ Batch, data = data, subset = Batch %in% c("Con", "Cis"))
test_LRY15_Cis <- wilcox.test(intrinsic_apoptosis_by_P531 ~ Batch, data = data, subset = Batch %in% c("LRY15", "Cis"))
test_Con_LRY15 <- wilcox.test(intrinsic_apoptosis_by_P531 ~ Batch, data = data, subset = Batch %in% c("LRY15", "Con"))
aggregate(intrinsic_apoptosis_by_P531 ~ Batch, data = data, mean)
ggplot(data, aes(Batch, intrinsic_apoptosis_by_P531, color = Batch))+geom_boxplot()+theme_bw()+RotatedAxis() + geom_text(aes(x = 1.5, y = 0.1, label = paste("p = ", format.pval(test_Con_Cis$p.value, digits = 3))), vjust = -0.5)+geom_text(aes(x = 2.5, y = 0.1, label = paste("p = ", format.pval(test_LRY15_Cis$p.value, digits = 3))), vjust = -0.5)+geom_text(aes(x = 2, y = 0.18, label = paste("p = ", format.pval(test_Con_LRY15$p.value, digits = 3))), vjust = -0.5) + scale_color_manual(values = c("#D55640","#6CB8D2","#E69F84"))+labs(x = "Oocyte", y = "intrinsic_apoptosis_by_P53")+theme(legend.background = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylim(-0.15, 0.3)
#Spatial gene set score 
Con_1_oocyte <- subset(oocyte, subset = sample == "Con_1")
data <- FetchData(Con1_oocyte, vars = c("x", "y", "intrinsic_apoptosis_by_P531"))
p_intrintic_apoptosis_by_P53 <- ggplot() +
  geom_point(
    data = Con_1_oocyte@meta.data[Con_1_oocyte@meta.data$sample == "Con_1", ],
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
      color = intrinsic_apoptosis_by_P531
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
    breaks = c(-0.2, 0, 0.2),
    limits = c(-0.2, 0.2),
    name = "Score"
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right"
  )
ggsave(
  filename = "Con_1_apoptosis_score_spatial.pdf",
  plot = p_intrintic_apoptosis_by_P53,
  width = 5,
  height = 4,
  units = "in"
)
#Cellchat
Con <- subset(Bin50, subset = Batch == "Con")
Cis <- subset(Bin50, subset = Batch == "Cis")
LRY15 <- subset(Bin50, subset = Batch == "LRY15")
#Con
data.input <- GetAssayData(Con, slot = "data")
meta <- data.frame(labels = Idents(seurat_object), row.names = colnames(data.input))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_Con.rds")
#Cis
data.input <- GetAssayData(Cis, slot = "data")
meta <- data.frame(labels = Idents(seurat_object), row.names = colnames(data.input))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_Cis.rds")
#LRY15
data.input <- GetAssayData(Cis, slot = "data")
meta <- data.frame(labels = Idents(seurat_object), row.names = colnames(data.input))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_LRY15.rds")
#Multiple_merge
cellchat_Con <- readRDS("./cellchat_Con.rds")
cellchat_Cis <- readRDS("./cellchat_Cis.rds")
cellchat_LRY15 <- readRDS("./cellchat_LRY15.rds")
object.list <- list(
  Con = cellchat_Con,
  Cis = cellchat_Cis,
  LRY15 = cellchat_LRY15
)
cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)
gg1 <- netVisual_bubble(cellchat_merged, sources.use = 1, targets.use = c(1:9),  comparison = c(2, 3), max.dataset = 3, title.name = "Increased signaling in LRY15", angle.x = 45, remove.isolate = T)
ggsave(
  filename = "bubble_Increased_signaling_in_LRY15.pdf",
  plot = gg1,
  width = 8,
  height = 6,
  units = "in"
)
pathways.show <- c("IGF")
weight.max <- getMaxWeight(
  object.list,
  slot.name = c("netP"),
  attribute = pathways.show
)
pdf(
  file = "CellChat_IGF_circle_Con_Cis_LRY15.pdf",
  width = 12,
  height = 4
)
par(mfrow = c(1, 3), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
dev.off()
#Con_1(Cis_2 and LRY15_6)
pdf(".Con_1_igf1_igf1r.pdf", width = 50, height = 10)
FeaturePlot(Con_1,features = c("Igf1", "Igf1r"),cols= c("grey", "red", "blue", "green"), blend = T, reduction = "spatial")
dev.off()
Bin50@meta.data$color_by_cluster <- ifelse(
  Bin50@meta.data$celltype == "Oocyte", "#f8d90d",
  ifelse(Bin50_merge@meta.data$cell_type == "Granulosa cell", "#00c6ff", "#E0EAFC")
)
Oocyte_granulosa_cell_Con_1spatial <- ggplot(
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
  filename = "./Oocyte_granulosa_cell_Con_1spatial.pdf",
  plot = Oocyte_granulosa_cell_Con_1spatial,
  width = 5,
  height = 4,
  units = "in"
)