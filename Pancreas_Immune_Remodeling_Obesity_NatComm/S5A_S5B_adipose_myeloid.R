# Pancreatic macrophage gene programs scored on the Massier adipose myeloid atlas
# (clusters myC0-myC15; Massier et al. 2023, GEO GSE225700).
# S5A: DotPlot of program genes. S5B: module-score heatmap.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

DATA_PATH  <- "path/to/object.rds"     # set to your local file
OUTPUT_DIR <- "path/to/output"         # set to your output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Load Massier myeloid object, RNA assay
my_seurat <- readRDS(DATA_PATH)
DefaultAssay(my_seurat) <- "RNA"

# Relabel clusters 0-15 to Massier myeloid groups
cluster_to_group <- c(
  "0"  = "M2 (myC0)",
  "1"  = "M2 (myC01)",
  "2"  = "LAM (myC02)",
  "3"  = "DC2 (myC03)",
  "4"  = "M2 (myC04)",
  "5"  = "non-classical Mo (myC05)",
  "6"  = "M1/M2-like (myC06)",
  "7"  = "M2 (myC07)",
  "8"  = "M2 (myC08)",
  "9"  = "M2 (myC09)",
  "10" = "MMe (myC10)",
  "11" = "M2 (myC11)",
  "12" = "M2 (myC12)",
  "13" = "M2 (myC13)",
  "14" = "classical Mo (myC14)",
  "15" = "Mox (myC15)"
)
Massier_group <- cluster_to_group[as.character(Idents(my_seurat))]
names(Massier_group) <- colnames(my_seurat)
my_seurat$Massier_group <- Massier_group

# S5A: DotPlot of macrophage program genes across myeloid clusters.
# Matches the published S5A panel as-is.
my_seurat$Massier_group <- factor(my_seurat$Massier_group, levels = c(
  "M2 (myC0)","M2 (myC01)","M2 (myC04)","M2 (myC13)","M2 (myC09)","M2 (myC08)",
  "M2 (myC12)","M2 (myC11)","M2 (myC07)","MMe (myC10)","LAM (myC02)","Mox (myC15)",
  "M1/M2-like (myC06)","DC2 (myC03)","classical Mo (myC14)","non-classical Mo (myC05)"
))

genes_paper <- c(
  "NUSAP1","HIST1H1B","CENPF","PCNA","MCM6",
  "CD58","YWHAH","HTRA1","BHLHE41","C3","CX3CR1",
  "PDPN","OLR1","APOE","CYP27A1","FABP5","LPL","APOC1",
  "NCEH1","CTSD","CD68","ACP5","PSAP","ASAH1","CD63","CTSA","MARCO",
  "LGALS3","LYZ","CD36","EREG","FCN1","VCAN","S100A12",
  "CD44","SELL","CCR2","CCR1","GAS6","AXL","MERTK",
  "CD209","CD163","MRC1","MSR1","MAF","PLTP","SELENOP","LYVE1","STAB1","FOLR2"
)
genes_paper <- genes_paper[genes_paper %in% rownames(my_seurat)]

# Summarized DotPlot data
p_tmp <- DotPlot(my_seurat, features = genes_paper, group.by = "Massier_group", dot.scale = 5)
df <- p_tmp$data
df$id <- factor(df$id, levels = levels(my_seurat$Massier_group))

# ggplot: shape-21 filled circles, blue-white-red scaled expression
p_dot_paper_genes <- ggplot(df, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Scaled\nexpression") +
  scale_size(range = c(0, 6), name = "Fraction\nin group (%)") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(face = "bold"),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    plot.margin     = margin(5, 5, 5, 5)
  )
ggsave(file.path(OUTPUT_DIR, "DotPlot_paper_genes_red_blue_blackframe.pdf"), plot = p_dot_paper_genes, width = 7, height = 9, dpi = 300)

# S5B: module-score heatmap.
# Code produces 10 programs x 16 myC columns; published S5B collapses the 16 clusters
# into 3 super-groups (Fetal-like / LAM / Monocytic); the 3-column version has no code.
gene_groups <- list(
  Fetal_derived_macrophage = c("FOLR2","STAB1","LYVE1","SELENOP","PLTP","MAF","MSR1","MRC1","CD163","CD209"),
  Endocytosis              = c("MERTK","AXL"),
  Efferocytosis            = c("MERTK","AXL","GAS6"),
  Monocyte_Migration       = c("CCR1","CCR2","SELL","CD44"),
  Monocyte_Markers         = c("S100A12","VCAN","FCN1","EREG","CD36","LYZ","LGALS3"),
  Lysosome                 = c("MARCO","CTSA","CD63","ASAH1","PSAP","ACP5","CD68","CTSD","NCEH1"),
  Lipid_Metabolism         = c("APOC1","LPL","FABP5","CYP27A1","APOE","OLR1","PDPN"),
  Pro_inflammatory_tissue_macrophages = c("CX3CR1","C3","BHLHE41","HTRA1","YWHAH","CD58"),
  Late_Mitosis             = c("MCM6","PCNA","CENPF"),
  Early_Mitosis            = c("HIST1H1B","NUSAP1")
)
gene_groups_filt <- lapply(gene_groups, function(g) g[g %in% rownames(my_seurat)])

my_seurat <- AddModuleScore(object = my_seurat, features = gene_groups_filt,
                            name = names(gene_groups_filt), assay = "RNA", nbin = 24, ctrl = 100)

score_cols <- c("Fetal_derived_macrophage1","Endocytosis2","Efferocytosis3","Monocyte_Migration4",
                "Monocyte_Markers5","Lysosome6","Lipid_Metabolism7",
                "Pro_inflammatory_tissue_macrophages8","Late_Mitosis9","Early_Mitosis10")

# Mean score per cluster, long format
scores_by_cluster <- my_seurat@meta.data %>%
  group_by(Massier_group) %>%
  summarise(across(all_of(score_cols), mean, .names = "{.col}")) %>%
  arrange(Massier_group)
scores_long <- scores_by_cluster %>%
  pivot_longer(cols = all_of(score_cols), names_to = "Module", values_to = "Score")
scores_long$Module <- gsub("[0-9]+$", "", scores_long$Module)

scores_long$Massier_group <- factor(scores_long$Massier_group, levels = c(
  "M2 (myC0)","M2 (myC01)","M2 (myC04)","M2 (myC13)","M2 (myC09)","M2 (myC08)",
  "M2 (myC12)","M2 (myC11)","M2 (myC07)","MMe (myC10)","LAM (myC02)","Mox (myC15)",
  "M1/M2-like (myC06)","DC2 (myC03)","classical Mo (myC14)","non-classical Mo (myC05)"
))
scores_long$Module <- factor(scores_long$Module, levels = c(
  "Early_Mitosis","Late_Mitosis","Pro_inflammatory_tissue_macrophages","Lipid_Metabolism",
  "Lysosome","Monocyte_Markers","Monocyte_Migration","Efferocytosis","Endocytosis",
  "Fetal_derived_macrophage"
))
scores_long <- dplyr::filter(scores_long, !is.na(Module))

p_module_heatmap <- ggplot(scores_long, aes(x = Massier_group, y = Module, fill = Score)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Module\nscore") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    plot.margin  = margin(5, 5, 5, 5)
  )
ggsave(file.path(OUTPUT_DIR, "ModuleScores_heatmap_CORRECTED.pdf"), plot = p_module_heatmap, width = 8, height = 5, dpi = 300)
