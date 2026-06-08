library(Seurat)
library(monocle)
library(Matrix)
library(Biobase)
library(ggplot2)
library(dplyr)
granulosa_cell <- readRDS("./Granulosa_cell.rds")
data <- as(as.matrix(Granuloda_cell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = data@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())	
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
express_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, express_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
saveRDS(cds,file="./monocle.rds")
Cell_type_pseudotime <- plot_cell_trajectory(monocds,color_by="seurat_clusters", size=1,show_backbone=TRUE)+ facet_wrap("~seurat_clusters", nrow = 1)+ scale_color_manual(values = colour)
ggsave(
  filename = "Cell_type_pseudotime.pdf",
  plot = Cell_type_pseudotime,
  width = 7,
  height = 2.5,
  units = "in"
)
#cell density pseudotime
pseudotime_data <- data.frame(
  Pseudotime = pData(monocds)$Pseudotime,
  Batch = pData(monocds)$Batch
)
pseudotime_data <- data.frame(
  Pseudotime = pData(monocds)$Pseudotime,
  Batch = pData(monocds)$Batch
)
p_pseudotime_density <- ggplot(
  pseudotime_data,
  aes(
    x = Pseudotime,
    fill = Batch
  )
) +
  geom_density(
    alpha = 0.9,
    color = "black",
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ Batch,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "CON" = "#D55640",
      "CIS" = "#6CB8D2",
      "LRY15" = "#E69F84"
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10)
  ) +
  scale_y_continuous(
    limits = c(0, 0.3),
    breaks = c(0, 0.1, 0.2, 0.3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Pseudotime",
    y = "Cell density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      color = "black",
      size = 14
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
ggsave(
  filename = "Cell_density_Pseudotime_density_CON_CIS_LRY15.pdf",
  plot = p_pseudotime_density,
  width = 7,
  height = 2.5,
  units = "in"
)
#Finding genes that changes as a function of pseudotime
Time_diff <- differentialGeneTest(
  cds[ordergene, ],
  cores = 1,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
Time_diff <- Time_diff[, c(5, 2, 3, 4, 1, 6, 7)]
write.csv(Time_diff,"Time_diff_all.csv",row.names = FALSE)
Time_genes <- Time_diff %>% arrange(qval) %>% slice_head(n = 100) %>% pull(gene_short_name) %>% as.character()
pdf(file = "Time_heatmapTop100.pdf", width = 5, height = 10)
p <- plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters = 4, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()
clusres <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.csv(clustering, "Time_clustering_all.csv", row.names = F)