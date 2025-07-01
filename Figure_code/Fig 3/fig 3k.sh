
args <- commandArgs(trailingOnly = TRUE)

effect_file <- args[1]
sign_file <- args[2]
output_file <- args[3]

data_effect <- read.table(effect_file, header = FALSE)
colnames(data_effect) <- c("Count_effect", "Value")
sorted_data_effect <- data_effect[order(data_effect$Value), ]
total_effect <- sum(sorted_data_effect$Count_effect)
sorted_data_effect$percentages_effect <- (sorted_data_effect$Count_effect / total_effect) * 100
selected_sorted_data_effect <- sorted_data_effect[, c("Value", "percentages_effect")]

data_sign <- read.table(sign_file, header = FALSE)
colnames(data_sign) <- c("Count_sign", "Value")
sorted_data_sign <- data_sign[order(data_sign$Value), ]
total_sign <- sum(sorted_data_sign$Count_sign)
sorted_data_sign$percentages_sign <- (sorted_data_sign$Count_sign / total_sign) * 100
selected_sorted_data_sign <- sorted_data_sign[, c("Value", "percentages_sign")]

merged_data <- merge(selected_sorted_data_sign, selected_sorted_data_effect, by = "Value", all = TRUE)
merged_data[is.na(merged_data)] <- 0

pdf(file = output_file, width = 8, height = 6)

ylim_fixed <- c(-20, 20)

barplot(
  merged_data$percentages_sign,
  names.arg = merged_data$Value,
  col = "orange",
  ylim = ylim_fixed,
  main = "Percentage Distribution (Up and Down)",
  xlab = "Value",
  ylab = "Percentage (%)",
  border = "black"
)

barplot(
  -merged_data$percentages_effect,
  col = "skyblue",
  add = TRUE,
  border = "black"
)

dev.off()
