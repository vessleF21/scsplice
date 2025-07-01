

library(ggplot2)

panel3 <- ggplot(data, aes(x = perc_covered_bases, fill = factor(bin_numReads))) + 
  theme_minimal() + 
  theme_classic() + 
  geom_density(alpha = .6) +
  scale_fill_brewer(palette = 'Set2') +
  xlab('Percentage of Covered Bases') + 
  ylab('Density') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_vline(xintercept = median(data$perc_covered_bases), linetype = 2, color = 'red')
ggsave('exon_con.pdf', panel3, width = 5, height = 3)