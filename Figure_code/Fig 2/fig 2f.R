f1 <- ggplot(CTLA4_data, aes(x=iid_label, y=celltype, fill=average_normalized_splicing)) + geom_tile() +
  scale_fill_viridis_c() + 
  theme(legend.position = 'top',
        axis.text.y = element_text(colour = celltype_pal[as.vector(ord_level2)]),
        axis.text.x = element_text(colour = ifelse(ord_level %in% iid_label[cluster.sig], 'red', 'grey'), angle = 45, hjust = 1))
ggsave('CTLA4.pdf', f1, width = 6, height = 6)