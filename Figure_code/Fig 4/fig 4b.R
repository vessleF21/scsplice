library(eulerr)

venn_data <- euler(c(
  "paSNP-gene pairs" = 7229,
  "pairs" = 518,
  "paSNP-gene pairs&pairs" = 692
))

plot(venn_data,
     fills = c("white", "gray50"),
     edges = "black",
     quantities = TRUE,
     labels = list(font = 2))
