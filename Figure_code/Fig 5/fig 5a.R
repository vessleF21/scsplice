rm(list=ls())
gc()

library(ggplot2)
library(reshape2)
library(dplyr)

ggf <- read.table("merged.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
ggf1 <- ggf[-30, ]

ggf1$Trait <- sub("^.*_(.*)$", "\\1", ggf1$Trait)

replacement_rules <- c(
  "Prostatecancer" = "PCa",
  "lungcarcinoma" = "LC",
  "Breastcancer" = "BrCa",
  "Basalcellcarcinoma" = "BCC",
  "hyperlipidemia" = "HPL",
  "Type2diabetes" = "T2D",
  "Hypertension" = "HTN",
  "Metabolicsyndrome" = "MetS"
)

ggf1 <- ggf1 %>%
  mutate(Trait = case_when(
    Trait %in% names(replacement_rules) ~ replacement_rules[Trait],
    TRUE ~ Trait
  ))

color_mapping <- c(
  "BNav" = "#01579B",
  "BMem" = "#A5EDFF",
  "CD4EM" = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157",
  "CD8GZMH" = "#FF5722",
  "CD8GZMK" = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58",
  "MonocM" = "#6A1B9A",
  "NoClaM" = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim" = "#D7CCC8"
)

colors_disease <- c(
  "Immune traits" = "#e78ac3",
  "Metabolic traits" = "#fc8d62",
  "Autoimmune diseases" = "#8da0cb",
  "Tumor diseases" = "#e5c494"
)

columns_to_divide <- c("BMem", "BNav", "CD4EM", "CD4Naive", "CD4Treg", 
                       "CD8GZMH", "CD8GZMK", "CD8Naive", "MAIT", "MonocM", 
                       "NKBright", "NKDim", "NoClaM")

for (col in columns_to_divide) {
  ggf1[[paste0(col, "_colocpr")]] <- ggf1[[col]] / ggf1$total_sum
}
for (col in columns_to_divide) {
  ggf1[[paste0(col, "_pr")]] <- ggf1[[col]] / ggf1$Total_GWAS_Loci
}

levels_order <- c("BAS","EOS","Hb","Ht","LYM","MCHC","MCH","MCV","MON","NEU","PLT","RBC","WBC","BMI","height","MetS","HTN","T2D","HPL","MS","SLE","RA","T1D","CD","UC","GD","BCC","BrCa","LC","PCa")
ggf1$Trait <- factor(ggf1$Trait, levels = rev(levels_order))

ggf1 <- ggf1 %>%
  mutate(
    disease_type = case_when(
      Trait %in% c("BAS","EOS","Hb","Ht","LYM","MCHC","MCH","MCV","MON","NEU","PLT","RBC","WBC") ~ "Immune traits",
      Trait %in% c("BMI","height","MetS","HTN","T2D","HPL") ~ "Metabolic traits",
      Trait %in% c("MS","SLE","RA","T1D","CD","UC","GD") ~ "Autoimmune diseases",
      Trait %in% c("BCC","BrCa","LC","PCa") ~ "Tumor diseases",
      TRUE ~ "Anthropometric traits"
    )
  )

ggff1 <- expand.grid(trait = ggf1$Trait, cell.type = columns_to_divide) %>%
  mutate(
    percent.of.colocalized.GWAS.loci = c(ggf1$BMem_pr,ggf1$BNav_pr,ggf1$CD4EM_pr,ggf1$CD4Naive_pr,ggf1$CD4Treg_pr,ggf1$CD8GZMH_pr,ggf1$CD8GZMK_pr,ggf1$CD8Naive_pr,ggf1$MAIT_pr,ggf1$MonocM_pr,ggf1$NKBright_pr,ggf1$NKDim_pr,ggf1$NoClaM_pr),
    number.of.colocalization = c(ggf1$BMem_colocpr,ggf1$BNav_colocpr,ggf1$CD4EM_colocpr,ggf1$CD4Naive_colocpr,ggf1$CD4Treg_colocpr,ggf1$CD8GZMH_colocpr,ggf1$CD8GZMK_colocpr,ggf1$CD8Naive_colocpr,ggf1$MAIT_colocpr,ggf1$MonocM_colocpr,ggf1$NKBright_pr,ggf1$NKDim_colocpr,ggf1$NoClaM_colocpr)
  )

ggff1$trait <- sub("^.*_(.*)$", "\\1", ggff1$trait)

ggff1$cell.type <- factor(ggff1$cell.type, levels = c("BNav","BMem","CD4Naive","CD4EM","CD4Treg",
                                                      "CD8Naive","CD8GZMH","CD8GZMK","MAIT", 
                                                      "NKDim", "NKBright","MonocM","NoClaM"))

ggff1$trait <- factor(ggff1$trait, levels = rev(levels_order))

ggplot(ggff1, aes(x = cell.type, y = trait)) +
  geom_tile(aes(fill = number.of.colocalization), color = "white", linewidth = 0.5) +
  geom_point(
    aes(size = percent.of.colocalized.GWAS.loci),
    shape = 21, color = "black", fill = "white", stroke = 0.5
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#7AA6DC", "#EFC000", "#CD534C"))(100),
    name = "Coloc Count"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    breaks = c(0.1, 0.2, 0.3),
    name = "Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey95", color = NA),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  ) +
  labs(x = NULL, y = NULL) +
  geom_hline(yintercept = seq(0.5, length(unique(ggff1$trait)) + 0.5, 1), 
             color = "gray80", linewidth = 0.2) +
  geom_vline(xintercept = seq(0.5, length(unique(ggff1$cell.type)) + 0.5, 1), 
             color = "gray80", linewidth = 0.2)
ggsave("colocalization_heatmap.pdf", width = 10, height = 10, dpi = 300, device = cairo_pdf)

colors_upper <- c(
  "#01579B", "#A5EDFF", "#1B5E20", "#8BC34A", "#D4E157", "#FF5722", "#FFCDD2", "#F9A825", "#FFEE58", "#6A1B9A", "#E1BEE7", "#82581F","#D7CCC8"
)

colors_disease <- c(
  "Immune traits" = "#8BC34A",
  "Metabolic traits" = "#fc8d62",
  "Autoimmune diseases" = "#8da0cb",
  "Tumor diseases" = "#e5c494"
)

cell_types <- c("BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
                "CD8Naive","CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "NKBright","MonocM","NoClaM")

column_sums <- colSums(ggf1[, 2:14])
column_sums_df <- data.frame(Column = names(column_sums), Sum = column_sums, stringsAsFactors = FALSE)
desired_order <- c("BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg", "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "NKBright", "MonocM", "NoClaM")
upper_df <- column_sums_df[match(desired_order, column_sums_df$Column), ]
colnames(upper_df) <- c("celltype", "num_coloc")
upper_df$celltype <- factor(upper_df$celltype, levels = desired_order)

ggplot(upper_df, aes(x = celltype, y = num_coloc, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = setNames(colors_upper, cell_types)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Cell Type Colocalization Counts")
ggsave("up.pdf", width = 10, height = 10, dpi = 300, device = cairo_pdf)

left_df <- data.frame(Column = ggf1$Trait, Sum = ggf1$total_sum, disease_type = ggf1$disease_type, stringsAsFactors = FALSE)

ggplot(left_df, aes(x = -Sum, y = Column, fill = disease_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = colors_disease) +
  scale_x_continuous(
    labels = function(x) abs(x),
    name = "Proportion of Colocalized Loci"
  ) +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    axis.text.y.right = element_text(face = "bold", hjust = 0),
    axis.ticks.y = element_line(),
    legend.position = "right"
  ) +
  labs(x = "Proportion of Colocalized Loci", y = NULL)
ggsave("left.pdf", width = 10, height = 10, dpi = 300, device = cairo_pdf)

right_df <- data.frame(Column = ggf1$Trait, Sum = ggf1$Total_GWAS_Loci, disease_type = ggf1$disease_type, stringsAsFactors = FALSE)

ggplot(right_df, aes(x = Column, y = Sum, fill = disease_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = colors_disease) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(y = "Proportion of Colocalized Loci", x = NULL)
ggsave("right.pdf", width = 10, height = 10, dpi = 300, device = cairo_pdf)
