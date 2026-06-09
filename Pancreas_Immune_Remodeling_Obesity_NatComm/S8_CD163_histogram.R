# CD163 expression across macrophage subsets vs Acinar/Islet negative controls,
# with a 95% lower-limit (5th-percentile) threshold. IMC, SpatialExperiment.

library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

DATA_PATH  <- "path/to/object.rds"     # set to your local file
OUTPUT_DIR <- "path/to/output"         # set to your output directory

# Load SpatialExperiment object
spe <- readRDS(DATA_PATH)

# CD163 marker, raw counts assay
cd163_marker <- rownames(spe)[grepl("CD163", rownames(spe), ignore.case = TRUE)]
cd163_expr   <- assay(spe, "counts")[cd163_marker, ]

# Macrophage subsets present in the data
MACROPHAGE_POPULATIONS <- c(
    "12_CD11cn_FOLR2_CD209_POS",
    "13_CD11cp_FOLR2_CD209_POS",
    "14_CD11cp_FOLR2_CD209_NEG",
    "15_CD11cn_FOLR2_CD209_NEG",
    "11_MPO_macrophages"
)
existing_macs <- MACROPHAGE_POPULATIONS[MACROPHAGE_POPULATIONS %in% unique(spe$celltype)]

# Per-population CD163
mac_cd163_data <- data.frame()
for (pop in existing_macs) {
    mac_cd163_data <- rbind(mac_cd163_data, data.frame(
        CD163 = cd163_expr[spe$celltype == pop], Population = pop, Category = "Macrophage"))
}

# All macrophages combined
all_macs   <- mac_cd163_data %>% filter(Category == "Macrophage") %>% mutate(Group = "All Macrophages")
mac_median <- median(all_macs$CD163, na.rm = TRUE)
n_macs     <- nrow(all_macs)

# Acinar (Gal-neg) and Islet negative controls
acinar_galneg <- data.frame(CD163 = cd163_expr[spe$celltype == "5_Galectin_NEG_acinar_cells"])
islets        <- data.frame(CD163 = cd163_expr[spe$celltype == "10_Islet_cells"])
acinar_galneg_median <- median(acinar_galneg$CD163, na.rm = TRUE)
islet_median         <- median(islets$CD163, na.rm = TRUE)
n_acinar_galneg <- nrow(acinar_galneg)
n_islets        <- nrow(islets)

# Per-subset CD163 values and plotting labels
mac_subsets <- list()
for (pop in existing_macs) mac_subsets[[pop]] <- cd163_expr[spe$celltype == pop]
short_names <- c(
    "12_CD11cn_FOLR2_CD209_POS" = "CD11c- FOLR2+",
    "13_CD11cp_FOLR2_CD209_POS" = "CD11c+ FOLR2+",
    "14_CD11cp_FOLR2_CD209_NEG" = "CD11c+ FOLR2-",
    "15_CD11cn_FOLR2_CD209_NEG" = "CD11c- FOLR2-",
    "11_MPO_macrophages"        = "MPO+ Monocytic",
    "Acinar_Galneg"             = "Acinar Gal-neg",
    "Islet_cells"               = "Islet Cells"
)

# Per-cell source data (CD163, Population)
detailed_source_data <- data.frame(CD163 = all_macs$CD163, Population = "All_Macrophages")
for (pop in existing_macs) {
    short_name <- ifelse(pop %in% names(short_names), short_names[pop], pop)
    detailed_source_data <- rbind(detailed_source_data, data.frame(
        CD163 = mac_subsets[[pop]], Population = gsub(" ", "_", short_name)))
}
detailed_source_data <- rbind(detailed_source_data, data.frame(CD163 = acinar_galneg$CD163, Population = "Acinar_Galneg"))
detailed_source_data <- rbind(detailed_source_data, data.frame(CD163 = islets$CD163, Population = "Islet_Cells"))

# 5th-percentile threshold (95% lower limit) over all macrophages
mac_95_lower     <- quantile(all_macs$CD163, probs = 0.05, na.rm = TRUE)
pct_acinar_below <- round(sum(acinar_galneg$CD163 < mac_95_lower) / n_acinar_galneg * 100, 1)
pct_islet_below  <- round(sum(islets$CD163 < mac_95_lower) / n_islets * 100, 1)

n_bins <- 100
threshold_plot_list <- list()

# All macrophages: 5th-percentile (blue), median (red)
threshold_plot_list[[1]] <- ggplot(data.frame(x = all_macs$CD163), aes(x = x)) +
    geom_histogram(bins = n_bins, fill = "#2ca02c", color = NA, alpha = 0.8, boundary = 0) +
    geom_vline(xintercept = mac_95_lower, linetype = "dashed", color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = mac_median, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = mac_95_lower + 3, y = Inf, label = paste0("5th %ile = ", round(mac_95_lower, 1)),
             vjust = 1.3, hjust = 0, color = "blue", size = 3.5, fontface = "bold") +
    annotate("text", x = mac_median + 3, y = Inf, label = paste0("Med = ", round(mac_median, 1)),
             vjust = 3, hjust = 0, color = "red", size = 3, fontface = "bold") +
    scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
    labs(title = paste0("All Macrophages (n=", format(n_macs, big.mark = ","), ")"), x = NULL, y = "Count") +
    theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 11))

# Per macrophage subset: % below the macrophage 5th percentile
for (i in seq_along(existing_macs)) {
    pop <- existing_macs[i]; vals <- mac_subsets[[pop]]
    med <- median(vals, na.rm = TRUE)
    pct_below  <- round(sum(vals < mac_95_lower) / length(vals) * 100, 1)
    short_name <- ifelse(pop %in% names(short_names), short_names[pop], pop)
    threshold_plot_list[[i + 1]] <- ggplot(data.frame(x = vals), aes(x = x)) +
        geom_histogram(bins = n_bins, fill = "#66c2a5", color = NA, alpha = 0.8, boundary = 0) +
        geom_vline(xintercept = mac_95_lower, linetype = "dashed", color = "blue", linewidth = 1.2) +
        geom_vline(xintercept = med, linetype = "dashed", color = "red", linewidth = 1) +
        annotate("text", x = mac_95_lower + 3, y = Inf, label = paste0(pct_below, "% below"),
                 vjust = 1.3, hjust = 0, color = "blue", size = 3, fontface = "bold") +
        annotate("text", x = med + 3, y = Inf, label = paste0("Med = ", round(med, 1)),
                 vjust = 3, hjust = 0, color = "red", size = 3, fontface = "bold") +
        scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
        labs(title = paste0(short_name, " (n=", format(length(vals), big.mark = ","), ")"), x = NULL, y = "Count") +
        theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10))
}

# Acinar Gal-neg control
threshold_plot_list[[length(threshold_plot_list) + 1]] <- ggplot(data.frame(x = acinar_galneg$CD163), aes(x = x)) +
    geom_histogram(bins = n_bins, fill = "#bdbdbd", color = NA, alpha = 0.8, boundary = 0) +
    geom_vline(xintercept = mac_95_lower, linetype = "dashed", color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = acinar_galneg_median, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = mac_95_lower + 3, y = Inf, label = paste0(pct_acinar_below, "% below mac 5th %ile"),
             vjust = 1.3, hjust = 0, color = "blue", size = 3.5, fontface = "bold") +
    annotate("text", x = 50, y = Inf, label = paste0("Med = ", round(acinar_galneg_median, 1)),
             vjust = 3, hjust = 0, color = "red", size = 3, fontface = "bold") +
    scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
    labs(title = paste0("Acinar Gal-neg (n=", format(n_acinar_galneg, big.mark = ","), ")"), x = NULL, y = "Count") +
    theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10))

# Islet control
threshold_plot_list[[length(threshold_plot_list) + 1]] <- ggplot(data.frame(x = islets$CD163), aes(x = x)) +
    geom_histogram(bins = n_bins, fill = "#d9d9d9", color = NA, alpha = 0.8, boundary = 0) +
    geom_vline(xintercept = mac_95_lower, linetype = "dashed", color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = islet_median, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = mac_95_lower + 3, y = Inf, label = paste0(pct_islet_below, "% below mac 5th %ile"),
             vjust = 1.3, hjust = 0, color = "blue", size = 3.5, fontface = "bold") +
    annotate("text", x = 50, y = Inf, label = paste0("Med = ", round(islet_median, 1)),
             vjust = 3, hjust = 0, color = "red", size = 3, fontface = "bold") +
    scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
    labs(title = paste0("Islet Cells (n=", format(n_islets, big.mark = ","), ")"), x = "CD163 Mean Intensity", y = "Count") +
    theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10))

# Stack panels
p_threshold <- wrap_plots(threshold_plot_list, ncol = 1) +
    plot_annotation(
        title = "CD163 Expression with 95% Lower Limit Threshold",
        subtitle = paste0("Blue line = macrophage 5th percentile (", round(mac_95_lower, 1), ") | Red line = median"),
        theme = theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 11)))

# Write figure
ggsave(file.path(OUTPUT_DIR, "Figure_CD163_95_LOWER_LIMIT_threshold.pdf"), p_threshold, width = 10, height = 14)

# Per-population threshold stats
threshold_stats <- data.frame(
    Population = c("All_Macrophages", existing_macs, "Acinar_Galneg", "Islet_cells"),
    n = c(n_macs, sapply(existing_macs, function(p) length(mac_subsets[[p]])), n_acinar_galneg, n_islets),
    median = c(mac_median, sapply(existing_macs, function(p) median(mac_subsets[[p]], na.rm = TRUE)), acinar_galneg_median, islet_median),
    mac_5th_percentile = mac_95_lower,
    pct_below_5th = c(5, sapply(existing_macs, function(p) round(sum(mac_subsets[[p]] < mac_95_lower) / length(mac_subsets[[p]]) * 100, 1)), pct_acinar_below, pct_islet_below)
)
write.csv(threshold_stats, file.path(OUTPUT_DIR, "Figure_CD163_95_LOWER_LIMIT_stats.csv"), row.names = FALSE)
write.csv(detailed_source_data, file.path(OUTPUT_DIR, "Figure_CD163_95_LOWER_LIMIT_SOURCE_DATA.csv"), row.names = FALSE)
