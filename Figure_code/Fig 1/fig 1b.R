library(tidyverse)
library(RColorBrewer)
rm(list=ls())
display.brewer.all()
c1<-brewer.pal(n=6,name="Set1")
pca_header = c('fid', 'iid', paste0('PC', 1:20))
pca_file = './merged.pca.eigenvec'
pca_result = read.table(pca_file, sep = '', header = F, stringsAsFactors = F, col.names = pca_header)
population = read.table('./sample-info.txt1', sep = '', header = F, stringsAsFactors = F, col.names = c('fid', 'iid', 'population'))
merged_df = left_join(pca_result, population, by=c('fid'='fid', 'iid'='iid'))
table(merged_df$population)
c_level = c('OneK1K', 'CEU', 'CHB', 'CHS', 'JPT','YRI')
eas_token = 1:6
# eas_token = c(2,3,4,5, 6,7)
# eas_token = c(1,2,4,5)
#eas_token = c(1,2)
c1 = c1[eas_token]
c_level = c_level[eas_token]
cols_t1 = c1[factor(merged_df$population, levels = c_level)]
pdf(paste0(pca_file,'.pdf') , family = 'serif')
# par(mfrow=c(2,2))
par(xpd = T, mar = par()$mar + c(0,0,0,3))
plot(PC2 ~ PC1, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC3 ~ PC2, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC4 ~ PC3, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC5 ~ PC4, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC6 ~ PC5, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC7 ~ PC6, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC8 ~ PC7, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC9 ~ PC8, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
plot(PC10 ~ PC9, merged_df, pch=21, bg=cols_t1, cex=0.7)
gap = (par('usr')[2] - par('usr')[1]) / 100
a = par('usr')[2] + gap
b = par('usr')[4] - (par('usr')[4] - par('usr')[3])/3
legend(a, b, legend=c_level, pt.bg= c1, pch=21, xpd=T)
par(mar = c(5,4,4,2) + 0.1)
dev.off()
