# Per-cell GZMB (Granzyme B) raw-intensity density for macrophages vs T cells,
# with an absolute positivity threshold of 2 (GZMB+ = raw intensity > 2). IMC, SpatialExperiment.

library(SpatialExperiment)
library(ggplot2)
library(patchwork)

DATA_PATH  <- "path/to/object.rds"     # set to your local file
OUTPUT_DIR <- "path/to/output"         # set to your output directory

# Load SpatialExperiment object
spe <- readRDS(DATA_PATH)

# Raw GZMB intensity per cell. Verify "Granzyme B Cell Intensity" matches a row of the object.
gzmb_raw <- assay(spe, "counts")["Granzyme B Cell Intensity", ]

# Macrophage and T-cell GZMB distributions
macs_data  <- data.frame(GZMB = gzmb_raw[spe$celltype == "12_CD11cn_FOLR2_CD209_POS"], CellType = "Macrophages")
tcell_data <- data.frame(GZMB = gzmb_raw[spe$celltype %in% c("8_CD4_T_cells", "9_CD8_T_cells")], CellType = "T_cells")

# Percent positive at threshold 2
macs_pct  <- round(mean(macs_data$GZMB > 2) * 100, 1)
tcell_pct <- round(mean(tcell_data$GZMB > 2) * 100, 1)

# Per-cell source data and summary stats
source_data <- rbind(macs_data, tcell_data)
summary_stats <- data.frame(
    CellType = c("Macrophages", "T_cells"),
    n_cells = c(nrow(macs_data), nrow(tcell_data)),
    mean_GZMB = c(mean(macs_data$GZMB), mean(tcell_data$GZMB)),
    median_GZMB = c(median(macs_data$GZMB), median(tcell_data$GZMB)),
    sd_GZMB = c(sd(macs_data$GZMB), sd(tcell_data$GZMB)),
    n_GZMB_pos = c(sum(macs_data$GZMB > 2), sum(tcell_data$GZMB > 2)),
    pct_GZMB_pos = c(macs_pct, tcell_pct),
    threshold = 2
)
write.csv(source_data, file.path(OUTPUT_DIR, "Figure_GZMB_Macs_vs_Tcells_source_data.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(OUTPUT_DIR, "Figure_GZMB_Macs_vs_Tcells_summary_stats.csv"), row.names = FALSE)

# Macrophage density, threshold line at 2
p_macs <- ggplot(macs_data, aes(x = GZMB)) +
    geom_density(fill = "#FF7F00", color = "black", alpha = 0.7, linewidth = 0.5) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 2.2, y = Inf, label = paste0("Threshold = 2\n", macs_pct, "% GZMB+"),
             vjust = 1.5, hjust = 0, color = "red", size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(0, 10), breaks = 0:10) +
    labs(title = paste0("Macrophages (n=", format(nrow(macs_data), big.mark = ","), ")"), x = NULL, y = "Density") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12))

# T-cell density, threshold line at 2
p_tcell <- ggplot(tcell_data, aes(x = GZMB)) +
    geom_density(fill = "#E41A1C", color = "black", alpha = 0.7, linewidth = 0.5) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 2.2, y = Inf, label = paste0("Threshold = 2\n", tcell_pct, "% GZMB+"),
             vjust = 1.5, hjust = 0, color = "red", size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(0, 10), breaks = 0:10) +
    labs(title = paste0("T Cells (n=", format(nrow(tcell_data), big.mark = ","), ")"), x = "GZMB Expression (raw counts)", y = "Density") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12))

# Stack and write
p_stacked <- p_macs / p_tcell +
    plot_annotation(title = "GZMB Distribution - Raw Counts",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(OUTPUT_DIR, "Figure_GZMB_Macs_vs_Tcells.pdf"), p_stacked, width = 10, height = 8)
