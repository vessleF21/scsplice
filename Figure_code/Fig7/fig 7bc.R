
rm(list = ls())
gc()

library(ggplot2)
library(ggpmisc)
library(dplyr)

COLORS_CELL <- c(
  "BNav"     = "#01579B",
  "BMem"     = "#A5EDFF",
  "CD4EM"    = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg"  = "#D4E157",
  "CD8GZMH"  = "#FF5722",
  "CD8GZMK"  = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT"     = "#FFEE58",
  "MonocM"   = "#6A1B9A",
  "NoClaM"   = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim"    = "#D7CCC8"
)

cell_summary <- read.table(
  "cell_type_summary_1.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""
)

plot_2 <- read.table(
  "plot_2.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""
)

df <- merge(cell_summary, plot_2, by = "cell_type")
colnames(df) <- c("cell_type", "n_donors", "mean_cells", "median_cells",
                  "total_cells", "mean_libsize", "median_libsize", "total_libsize",
                  "sgene", "mean_slope", "median_slope", "color",
                  "n_cis", "n_trans", "n_donors_y")

df$color <- COLORS_CELL[df$cell_type]

pearson_summary <- function(x, y, label) {
  ct  <- cor.test(x, y, method = "pearson")
  sig <- cut(ct$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
             labels = c("***", "**", "*", "ns"))
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("r = %.4f,  p = %s,  95%% CI: [%.4f, %.4f],  %s\n",
              ct$estimate, format.pval(ct$p.value, digits = 3),
              ct$conf.int[1], ct$conf.int[2], sig))
  ct
}

make_label <- function(ct) {
  paste0("r = ", round(ct$estimate, 3), "\n",
         "p = ", format.pval(ct$p.value, digits = 3), "\n",
         "95% CI: [", round(ct$conf.int[1], 3), ", ",
         round(ct$conf.int[2], 3), "]")
}

theme_scatter <- function() {
  theme_classic() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title       = element_text(face = "bold", size = 16, color = "black"),
      axis.text        = element_text(size = 14, color = "black"),
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6)
    )
}

make_scatter <- function(data, x_var, y_var, x_lab, cor_result) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = cell_type), size = 5, shape = 21,
               fill = "white", stroke = 2) +
    scale_color_manual(values = data$color) +
    scale_y_log10() +
    stat_poly_line(formula = y ~ x, color = "darkblue", linewidth = 1.5) +
    stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label)),
                 size = 4, label.x = 0.05, label.y = 0.95, color = "darkred") +
    annotate("text", x = Inf, y = Inf, label = make_label(cor_result),
             hjust = 1.1, vjust = 1.5, size = 5,
             color = "darkgreen", fontface = "bold") +
    theme_scatter() +
    xlab(x_lab) +
    ylab("Number of trans-sGenes") +
    guides(color = "none")
}

ct1 <- pearson_summary(df$n_donors, df$n_trans, "n_donors vs n_trans")
p1  <- make_scatter(df, "n_donors", "n_trans", "Number of donors", ct1)
print(p1)
ggsave("transgene_vs_donor.pdf", p1, height = 7, width = 7)

ct2 <- pearson_summary(df$mean_libsize, df$n_trans, "mean_libsize vs n_trans")
p2  <- make_scatter(df, "mean_libsize", "n_trans", "Junction read counts", ct2)
print(p2)
ggsave("transgene_vs_mean_libsize.pdf", p2, height = 7, width = 7)

