
rm(list = ls())
gc()

library(readxl)
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)
library(ggplot2)


colors2 <- c(
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

label_dict <- list(
  "BMem" = "BMem",
  "BNav" = "BNav",
  "CD4EM" = "CD4EM",
  "CD4Naive" = "CD4Naive",
  "CD4Treg" = "CD4Treg",
  "CD8GZMH" = "CD8GZMH",
  "CD8GZMK" = "CD8GZMK",
  "CD8Naive" = "CD8Naive",
  "MAIT" = "MAIT",
  "NKBright" = "NKBright",
  "NKDim" = "NKDim",
  "MonocM" = "MonocM",
  "NoClaM" = "NoClaM"
)

data<-read.csv("AS",sep='\t')

dat.filter<-data
for (i in seq_along(label_dict)) {
  dat.filter$Celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat.filter$Celltype)
}
summary <- dat.filter %>% group_by(Celltype) %>% summarise(median_n = median(num))
dat.filter$Celltype <- factor(dat.filter$Celltype, levels = summary$Celltype[order(summary$median_n, decreasing = T)])

pdf("as.pdf")

ggplot(dat.filter, aes(x=Celltype, y=num, fill=Celltype)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values=colors2) +
  theme_minimal() +
  theme_classic() +
  xlab("Cell type") + 
  ylab("Number of AS genes") +
  theme(legend.position = "none") +
  theme(    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")) +
  scale_y_continuous(breaks = c(1000, 2000, 3000, 4000, 5000, 6000))

dev.off()




