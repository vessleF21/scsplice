




rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(data.table)
library(scales)

WORK_DIR       <- "summary"
PLOT_DATA_PATH <- "plot_data_1.txt"

CELL_ORDER <- c(
  "BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "MonocM", "NoClaM"
)

TRAIT_ORDER <- c(
  "BAS", "EOS", "Hb", "Ht", "LYM", "MCHC", "MCH", "MCV", "MON", "NEU", "PLT", "RBC", "WBC",
  "BMI", "height", "MetS", "HTN", "T2D", "HPL",
  "MS", "SLE", "RA", "T1D", "CD", "UC", "GD",
  "BCC", "BrCa", "LC", "PCa"
)

TRAIT_RENAME <- c(
  "Prostatecancer"     = "PCa",
  "lungcarcinoma"      = "LC",
  "Breastcancer"       = "BrCa",
  "Basalcellcarcinoma" = "BCC",
  "hyperlipidemia"     = "HPL",
  "Type2diabetes"      = "T2D",
  "Hypertension"       = "HTN",
  "Metabolicsyndrome"  = "MetS"
)

DISEASE_TYPE_MAP <- list(
  "Immune traits"       = c("BAS","EOS","Hb","Ht","LYM","MCHC","MCH","MCV","MON","NEU","PLT","RBC","WBC"),
  "Metabolic traits"    = c("BMI","height","MetS","HTN","T2D","HPL"),
  "Autoimmune diseases" = c("MS","SLE","RA","T1D","CD","UC","GD"),
  "Tumor diseases"      = c("BCC","BrCa","LC","PCa")
)

COLORS_CELL <- c(
  "BNav"     = "#01579B",
  "BMem"     = "#A5E6FF",
  "CD4Naive" = "#8BC34A",
  "CD4EM"    = "#1B5E20",
  "CD4Treg"  = "#D4E157",
  "CD8Naive" = "#F9A825",
  "CD8GZMH"  = "#FF5722",
  "CD8GZMK"  = "#FFCDD2",
  "MAIT"     = "#FFEE58",
  "NKDim"    = "#D7CCC8",
  "NKBright" = "#82581F",
  "MonocM"   = "#6A1B9A",
  "NoClaM"   = "#E1BEE7"
)

COLORS_DISEASE <- c(
  "Immune traits"       = "#8BC34A",
  "Metabolic traits"    = "#fc8d62",
  "Autoimmune diseases" = "#8da0cb",
  "Tumor diseases"      = "#e5c494"
)

setwd(WORK_DIR)

raw <- read.table("merged_data_2.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
df  <- raw[-30, ]

df$Trait <- sub("^.*_(.*)$", "\\1", df$Trait)
df$Trait <- dplyr::recode(df$Trait, !!!TRAIT_RENAME)

df <- df %>%
  mutate(disease_type = case_when(
    Trait %in% DISEASE_TYPE_MAP[["Immune traits"]]       ~ "Immune traits",
    Trait %in% DISEASE_TYPE_MAP[["Metabolic traits"]]    ~ "Metabolic traits",
    Trait %in% DISEASE_TYPE_MAP[["Autoimmune diseases"]] ~ "Autoimmune diseases",
    Trait %in% DISEASE_TYPE_MAP[["Tumor diseases"]]      ~ "Tumor diseases",
    TRUE ~ "Other"
  ))

for (cell in CELL_ORDER) {
  df[[paste0(cell, "_colocpr")]] <- df[[cell]] / df$total_sum
  df[[paste0(cell, "_pr")]]      <- df[[cell]] / df$Total_GWAS_Loci
}

df$Trait <- factor(df$Trait, levels = rev(TRAIT_ORDER))

heatmap_df <- lapply(CELL_ORDER, function(cell) {
  data.frame(
    trait      = df$Trait,
    cell.type  = cell,
    pct_gwas   = df[[paste0(cell, "_pr")]],
    norm_coloc = df[[paste0(cell, "_colocpr")]],
    stringsAsFactors = FALSE
  )
}) |> do.call(what = rbind)

sqtl_ref <- read.table(PLOT_DATA_PATH, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "")

cell_sgene <- sqtl_ref %>%
  group_by(cell_type) %>%
  summarise(sgene = n(), .groups = "drop") %>%
  mutate(pro = sgene / 28) %>%
  rename(cell.type = cell_type)

heatmap_df <- heatmap_df %>%
  left_join(cell_sgene, by = "cell.type") %>%
  mutate(norm_coloc = norm_coloc / pro) %>%
  select(trait, cell.type, pct_gwas, norm_coloc)

heatmap_df$cell.type <- factor(heatmap_df$cell.type, levels = CELL_ORDER)
heatmap_df$trait     <- factor(heatmap_df$trait,     levels = rev(TRAIT_ORDER))

p_heatmap <- ggplot(heatmap_df, aes(x = cell.type, y = trait)) +
  geom_tile(aes(fill = norm_coloc), color = "white", linewidth = 0.5) +
  geom_point(aes(size = pct_gwas), shape = 21, color = "black", fill = "white", stroke = 0.5) +
  scale_fill_gradientn(
    colours = c("#BBDEFB", "#42A5F5", "#1976D2", "#0D47A1", "#002171"),
    limits  = c(0, 0.10), oob = squish,
    breaks  = seq(0, 0.10, 0.02),
    labels  = paste0(seq(0, 10, 2), "%"),
    name    = "Normalized\nColoc Count",
    guide   = guide_colorbar(barheight = unit(1.5, "cm"))
  ) +
  scale_size_continuous(
    range  = c(1, 8),
    breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.20),
    labels = paste0(c(0, 4, 8, 12, 16, 20), "%"),
    name   = "Proportion of colocalized GWAS loci"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_hline(yintercept = seq(0.5, length(unique(heatmap_df$trait))     + 0.5, 1), color = "gray80", linewidth = 0.2) +
  geom_vline(xintercept = seq(0.5, length(unique(heatmap_df$cell.type)) + 0.5, 1), color = "gray80", linewidth = 0.2) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y       = element_text(face = "bold", color = "black"),
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "right",
    legend.key.height = unit(1.5, "cm")
  ) +
  labs(x = NULL, y = NULL)

upper_df <- data.frame(
  celltype  = factor(CELL_ORDER, levels = CELL_ORDER),
  num_coloc = colSums(df[, CELL_ORDER])
)

p_upper <- ggplot(upper_df, aes(x = celltype, y = num_coloc, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = COLORS_CELL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.title      = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Cell Type Colocalization Counts")

left_df <- df %>%
  select(Trait, total_sum, disease_type) %>%
  mutate(Trait = factor(Trait, levels = rev(TRAIT_ORDER)))

p_left <- ggplot(left_df, aes(x = -total_sum, y = Trait, fill = disease_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = COLORS_DISEASE) +
  scale_x_continuous(labels = function(x) abs(x)) +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    axis.text.y.right = element_text(face = "bold", hjust = 0),
    axis.ticks.y      = element_line(),
    legend.position   = "right"
  ) +
  labs(x = "Number of Colocalized Loci", y = NULL)

right_df <- df %>%
  select(Trait, Total_GWAS_Loci, disease_type) %>%
  mutate(Trait = factor(Trait, levels = rev(TRAIT_ORDER)))

p_right <- ggplot(right_df, aes(x = Total_GWAS_Loci, y = Trait, fill = disease_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = COLORS_DISEASE) +
  theme_classic() +
  theme(
    axis.text.y     = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(x = "Total GWAS Loci", y = NULL)

print(p_heatmap)
print(p_upper)
print(p_left)
print(p_right)

ggsave("fig_heatmap.pdf", plot = p_heatmap, width = 10, height = 10)
ggsave("fig_upper.pdf",   plot = p_upper,   width = 10, height =  4)
ggsave("fig_left.pdf",    plot = p_left,    width =  5, height = 10)
ggsave("fig_right.pdf",   plot = p_right,   width =  5, height = 10)




