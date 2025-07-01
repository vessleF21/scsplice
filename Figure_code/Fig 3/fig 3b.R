

distance<-read.table("all",sep=" ",head=F)

frequency<-as.data.frame(matrix(0,110,2))

for(i in 1:nrow(distance)){
  if((distance[i,5]>=distance[i,3]) && (distance[i,5]<=distance[i,4])){
    pro <- (distance[i,5] - distance[i,3]) / (distance[i,4] - distance[i,3])
    print(pro)
    print(i)
    if(pro == 0){ pro <- "0.0" }
    
    pro_num <- as.numeric(substr(as.character(pro), 3, 3))
    
    if(!is.na(pro_num)){
      if(distance[i,2] == "+"){
        frequency[51 + pro_num, 1] <- frequency[51 + pro_num, 1] + 1
      }
      if(distance[i,2] == "-"){
        frequency[60 - pro_num, 1] <- frequency[60 - pro_num, 1] + 1
      }
    }
  }
  
  if(distance[i,5] < distance[i,3]){
    pro <- (distance[i,3] - distance[i,5]) / 500
    pro <- floor(pro)
    if(pro < 50){
      if(distance[i,2] == "+"){
        frequency[50 - pro, 1] <- frequency[50 - pro, 1] + 1
      }
      if(distance[i,2] == "-"){
        frequency[61 + pro, 1] <- frequency[61 + pro, 1] + 1
      }
    }
  }
  
  if(distance[i,5] > distance[i,4]){
    pro <- (distance[i,5] - distance[i,4]) / 500
    pro <- floor(pro)
    if(pro < 50){
      if(distance[i,2] == "+"){
        frequency[61 + pro, 1] <- frequency[61 + pro, 1] + 1
      }
      if(distance[i,2] == "-"){
        frequency[50 - pro, 1] <- frequency[50 - pro, 1] + 1
      }
    }
  }
}


library(ggplot2)
for(i in 1:50){
  frequency[51-i,2]<-paste(i,"*0.5kb upstream",sep="")
  frequency[60+i,2]<-paste(i,"*0.5kb downstream",sep="")
}
for(i in 1:10){
  frequency[50+i,2]<-paste(i," internal",sep="")
}
colnames(frequency)<-c("Number of sQTL","Position")

frequency$Position <- factor(frequency$Position, levels=frequency$Position)

pp <- ggplot(frequency, aes(x = Position, y = `Number of sQTL`)) +
  geom_bar(stat = "identity", fill = "#0066CC", color = NA, width = 0.6, alpha = 0.8) +  # Remove black borders
  geom_vline(xintercept = c(50, 61), linetype = "dashed", color = "red", size = 1.5) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", color = "#333333", hjust = 0.5, margin = margin(b = 20)),
    plot.subtitle = element_text(size = 16, color = "#666666", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold", color = "#1E1E1E"),
    axis.text = element_text(size = 14, color = "#1E1E1E"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "#333333"),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.7),
    legend.position = "none",
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"),
    plot.background = element_rect(fill = "white", color = NA),  # White background for the plot
    panel.background = element_blank()  # Remove background for the panel area
  ) +
  xlab("Position") +
  ylab("Number of sQTL") +
  ggtitle("Distribution of sQTL Across Position") +
  annotate("text", x = 50, y = 80, label = "Upstream region", angle = 0, size = 5, color = "#D32F2F", fontface = "italic") +
  annotate("text", x = 61, y = 80, label = "Downstream region", angle = 0, size = 5, color = "#D32F2F", fontface = "italic")

print(pp)
ggsave("position.pdf", pp, height = 10, width = 20)
