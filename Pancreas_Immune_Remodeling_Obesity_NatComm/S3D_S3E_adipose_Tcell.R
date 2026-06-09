# Pancreatic T/NK gene programs scored on the Massier adipose lymphocyte atlas
# (clusters lyC0-lyC10; Massier et al. 2023, GEO GSE225700).
# S3D: DotPlot of gene programs. S3E: module-score heatmap.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

DATA_PATH  <- "path/to/object.rds"     # set to your local file
OUTPUT_DIR <- "path/to/output"         # set to your output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Load Massier lymphoid object, RNA assay
ly_seurat <- readRDS(DATA_PATH)
DefaultAssay(ly_seurat) <- "RNA"

# Relabel clusters 0-10 to Massier lymphoid groups (Massier et al. 2023, Fig 3d/e)
cluster_to_group <- c(
  "0"  = "CD4+ Th1 (lyC0)",
  "1"  = "CD4+ TRM (lyC01)",
  "2"  = "CD8+ early-diff (lyC02)",
  "3"  = "CD8+ late-diff (lyC03)",
  "4"  = "CD4+ naive (lyC04)",
  "5"  = "NKT (lyC05)",
  "6"  = "NK CD16+ (lyC06)",
  "7"  = "lyC07",
  "8"  = "Treg (lyC08)",
  "9"  = "NK CD16- (lyC09)",
  "10" = "lyC10"
)
Massier_group <- cluster_to_group[as.character(Idents(ly_seurat))]
names(Massier_group) <- colnames(ly_seurat)
ly_seurat$Massier_group <- Massier_group
ly_seurat$Massier_group <- factor(ly_seurat$Massier_group, levels = c(
  "CD4+ Th1 (lyC0)", "CD4+ TRM (lyC01)", "CD4+ naive (lyC04)", "Treg (lyC08)",
  "CD8+ early-diff (lyC02)", "CD8+ late-diff (lyC03)", "NKT (lyC05)",
  "NK CD16+ (lyC06)", "NK CD16- (lyC09)", "lyC07", "lyC10"
))

# S3D: DotPlot of gene programs across lymphoid clusters.
# Code plots the full gene set (adds a 25-gene Obesity-IFN block, Fig 3E);
# published S3D shows the 3-program subset: Chemotaxis, Cytotoxic function, Negative regulation.
# CCL5/TCF7/TNFAIP3 UMAP feature plots on the same figure page have no saved code.
genes_trm_paper <- c(
  # Obesity-driven IFN/inflammatory response (Fig 3E)
  "NFKB2", "DDX58", "TRAF1", "KLF2", "G0S2", "S1PR4", "IFITM1",
  "HLA-DRA", "SOCS3", "VEGFA", "CSF3", "METRNL", "CMPK2", "TYROBP", "SOD2",
  "TIMP1", "ICAM1", "IL6", "CD74", "PLSCR1", "IFI30", "CD7",
  "MX2", "MX1", "STAT3",
  # Negative regulation
  "NR3C1", "BTG2", "NFKBIA", "ANXA1", "DUSP1", "BTG1", "ZFP36",
  "ZFP36L2", "TNFAIP3", "TSC22D3",
  # Cytotoxic function
  "NKG7", "IL2RG", "GZMH", "GZMB", "FASLG", "CTSW", "GZMA",
  "IL32", "GLUL", "PRF1", "CD2", "KLRC1", "ZNF683",
  # Chemotaxis
  "CCR9", "CXCR4", "GPR183", "CCR5", "CCR1", "CCR2",
  "CXCR3", "CXCR6", "GPR25"
)
genes_trm_paper <- genes_trm_paper[genes_trm_paper %in% rownames(ly_seurat)]

# Summarized DotPlot data
p_tmp <- DotPlot(ly_seurat, features = genes_trm_paper, group.by = "Massier_group", dot.scale = 5)
df <- p_tmp$data
df$id <- factor(df$id, levels = levels(ly_seurat$Massier_group))

# Source-data CSV
df_export <- df[, c("id", "features.plot", "avg.exp", "avg.exp.scaled", "pct.exp")]
colnames(df_export) <- c("Cluster", "Gene", "Avg_Expression", "Scaled_Expression", "Pct_Expressed")
write.csv(df_export, file = file.path(OUTPUT_DIR, "SourceData_DotPlot_TRM_genes.csv"), row.names = FALSE)

# ggplot: shape-21 filled circles, blue-white-red scaled expression
p_dot_pub <- ggplot(df, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Scaled\nexpression") +
  scale_size(range = c(0, 6), name = "Fraction\nin group (%)") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y     = element_text(face = "bold", size = 11),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  )
ggsave(file.path(OUTPUT_DIR, "DotPlot_TRM_genes_red_blue_blackframe.pdf"), plot = p_dot_pub, width = 8, height = 16, dpi = 300)

# S3E: module-score heatmap.
# Code scores 5 modules (adds Costimulatory and a 25-gene Obesity-IFN module);
# published S3E shows the 3-program subset: Cytotoxic function, Chemotaxis, Negative regulation.
gene_groups <- list(
  TRM_Cytotoxic_Effector = c("ZNF683","KLRC1","CD2","PRF1","GLUL","IL32",
                             "GZMA","CTSW","FASLG","GZMB","GZMH","IL2RG","NKG7"),
  TRM_Chemotaxis        = c("GPR25","CXCR6","CXCR3","CCR2","CCR1","CCR5",
                            "GPR183","CXCR4","CCR9"),
  TRM_Costimulatory     = c("CD2","CD28","ICOS"),
  TRM_Immunoreg_Brakes  = c("TSC22D3","TNFAIP3","ZFP36L2","ZFP36","BTG1",
                            "DUSP1","ANXA1","NFKBIA","BTG2","NR3C1"),
  TRM_Obesity_IFN       = c("STAT3","MX1","MX2","CD7","IFI30","PLSCR1","CD74","IL6",
                            "ICAM1","TIMP1","SOD2","TYROBP","CMPK2","METRNL","CSF3",
                            "VEGFA","SOCS3","HLA-DRA","IFITM1","S1PR4","G0S2","KLF2",
                            "TRAF1","DDX58","NFKB2")
)
gene_groups_filt <- lapply(gene_groups, function(g) g[g %in% rownames(ly_seurat)])
gene_groups_filt <- gene_groups_filt[sapply(gene_groups_filt, length) >= 2]

ly_seurat <- AddModuleScore(object = ly_seurat, features = gene_groups_filt,
                            name = names(gene_groups_filt), assay = "RNA", nbin = 24, ctrl = 100)
score_cols <- grep(paste0("^(", paste(names(gene_groups_filt), collapse = "|"), ")"),
                   colnames(ly_seurat@meta.data), value = TRUE)

# Mean score per cluster, long format
scores_by_cluster <- ly_seurat@meta.data %>%
  group_by(Massier_group) %>%
  summarise(n_cells = n(), across(all_of(score_cols), mean, .names = "{.col}")) %>%
  arrange(Massier_group)
scores_long <- scores_by_cluster %>% select(-n_cells) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "Module", values_to = "Score")
scores_long$Module <- gsub("[0-9]+$", "", scores_long$Module)

module_order <- rev(c("TRM_Cytotoxic_Effector","TRM_Chemotaxis","TRM_Costimulatory",
                      "TRM_Immunoreg_Brakes","TRM_Obesity_IFN"))
module_order <- module_order[module_order %in% unique(scores_long$Module)]
scores_long$Module <- factor(scores_long$Module, levels = module_order)

p_module_heatmap <- ggplot(scores_long, aes(x = Massier_group, y = Module, fill = Score)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Module\nscore") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y  = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(5, 5, 5, 5)
  ) +
  ggtitle("Pancreatic TRM Signatures on Adipose Lymphocyte Clusters")
ggsave(file.path(OUTPUT_DIR, "TRM_ModuleScores_heatmap.pdf"), plot = p_module_heatmap, width = 8, height = 4, dpi = 300)
