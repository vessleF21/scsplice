cc_trait<-c("ValetteK_CommunBiol_2021_asthma",
"LiuZ_NatGenet_2023_CD",
"AndlauerTF_SciAdv_2016_MS",
"Ishigaki_NatGenet_2022_RA",
"BenthamJ_NatGenet_2015_SLE",
"SakaueS_NatGenet_2021_SS",
"ChiouJ_Nature_2021_T1D",
"LiuZ_NatGenet_2023_UC",
"SakaueS_NatGenet_2021_GD",
"MichailidouK_Nature_2017_Breastcancer",
"WangA_NatGenet_2023_Prostatecancer",
"McKayJD_NatGenet_2017_lungcarcinoma",
"PurdueMP_NatGenet_2024_Kidneycancer",
"SeviiriM_NatCommun_2022_Basalcellcarcinoma",
"MaraTA_NatCommun_2018_endometrialcarcinoma",
"SeviiriM_NatCommun_2022_squamouscellcarcinoma",
"RamachandranD_HumMolGenet_2022_cervicalcancer",
"DarengEO_AmJHumGenet_2024_ovariancancer",
"SakaueS_NatGenet_2021_Type2diabetes",
"TrinderM_Atherosclerosis_2021_hyperlipidemia"
)

cc_trait_case<-c(56167,20873,4888,22350,5201,1599,18942,23252,4487,76192,122188,29266,29020,20791,8758,7402,363,15588,84224,17485)
cc_trait_control<-c(352255,346719,10395,74823,9066,658316,501638,352256,629598,63082,604640,56450,835670,286893,46126,286892,861,105724,583280,331737)

quant_trait <- c("ChenMH_CELL_2020_BAS",
"SakaueS_NatGenet_2021_BMI",
"ChenMH_CELL_2020_EOS",
"SakaueS_NatGenet_2021_height",
"ChenMH_CELL_2020_Ht",
"ChenMH_CELL_2020_LYM",
"ChenMH_CELL_2020_MCHC",
"ChenMH_CELL_2020_MCH",
"ChenMH_CELL_2020_MCV",
"ChenMH_CELL_2020_MON",
"ChenMH_CELL_2020_NEU",
"ChenMH_CELL_2020_PLT",
"ChenMH_CELL_2020_RBC",
"ChenMH_CELL_2020_WBC",
"ChenMH_CELL_2020_Hb",
"ParkS_NatGenet_2024_Metabolicsyndrome",
"CareyCE_NatHumBehav_2024_Hypertension"
)

quant_trait_num<-c(474001,523818,474237,525444,562259,524923,491553,486823,544127,521594,519288,542827,545203,562243,563946,1384348,338391)

cell_type<-c("BMem", "BNav", "CD4EM", "CD4Naive", "CD4Treg", "CD8GZMH", "CD8GZMK", "CD8Naive", "MAIT", "MonocM", "NKBright", "NKDim", "NoClaM")
celltype_sample<-c(865,873,980,975,729,939,924,759,352,603,346,976,435)

args<-(commandArgs(TRUE))
gwas_name<-args[1]
celltype<-args[2]
chr<-args[3]
intron<-args[4]

num_1<-which(cc_trait==gwas_name)
num_2<-which(quant_trait==gwas_name)
if(length(num_1)>0)trait_type<-"cc"
if(length(num_2)>0)trait_type<-"quant"
num_3<-which(cell_type==celltype)

library("coloc")
library(dplyr)
library(data.table)

GWAS_data <- paste0(celltype, "_", chr, "_nominals_1/", gwas_name, "/", intron)
print("start read in")
GWAS_pre<-fread(GWAS_data,sep="\t",header=TRUE)
GWAS_pre<-as.data.frame(GWAS_pre)
GWAS_pre<-GWAS_pre[,2:5]
colnames(GWAS_pre)<-c("variant_id","beta","se","pval_nominal")

MAF <- read.table(paste0(chr, ".txt"), header = TRUE)
colnames(MAF)<-c("variant_id","maf")
print("finished prepare")
test_rs<-as.data.frame(matrix(NA,1,16))
colnames(test_rs)<-c("sqtl_variant_id_most","intron_id","gene","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","sqtl_variant_P_most","sqtl_variant_id_merge","sqtl_variant_P_merge","GWAS_variant_id_most","GWAS_variant_P_most","GWAS_variant_id_merge","GWAS_variant_P_merge")


snp_1mb <- fread(paste0(celltype, "_", chr, "_nominals_1/", gwas_name, "/", intron), sep = " ", header = FALSE)
colnames(snp_1mb)<-c("gene","intron","variant_id","pval_nominal","beta")
snp_1mb<-as.data.frame(snp_1mb)
snp_1mb<-snp_1mb[order(as.numeric(snp_1mb$pval_nominal)),]

test_rs[1,1]<-snp_1mb[1,3]
test_rs[1,2]<-snp_1mb[1,2]
test_rs[1,3]<-snp_1mb[1,1]
test_rs[1,10]<-snp_1mb[1,4]

GWAS_pre_order <- GWAS_pre[order(as.numeric(GWAS_pre$pval_nominal)), ]
test_rs[1,13]<-GWAS_pre_order[1,1]
test_rs[1,14]<-GWAS_pre_order[1,4]

input<-merge(snp_1mb,MAF,by="variant_id")
input<-merge(input,GWAS_pre,by="variant_id",suffixes=c("_sqtl","_gwas"))
input<-input[order(as.numeric(input$pval_nominal_gwas)),]


input_ordersqtl<-input[order(as.numeric(input$pval_nominal_sqtl)),]
test_rs[1,11]<-input_ordersqtl[1,1]
test_rs[1,12]<-input_ordersqtl[1,4]

input_ordergwas<-input[order(as.numeric(input$pval_nominal_gwas)),]
test_rs[1,15]<-input_ordergwas[1,1]
test_rs[1,16]<-input_ordergwas[1,9]

input$varbeta<-(input$se)^2

if(nrow(input)>10){
        if(trait_type=="cc"){
                result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="cc", beta=input$beta_gwas,varbeta=input$varbeta,
                                                  s=cc_trait_case[num_1]/(cc_trait_case[num_1]+cc_trait_control[num_1]),
                                                  N=cc_trait_case[num_1]+cc_trait_control[num_1],snp = input$variant_id),
                                    dataset2=list(pvalues=input$pval_nominal_sqtl, type="quant", N=celltype_sample[num_3],snp=input$variant_id), MAF=input$maf)
                test_rs[1,4:9]<-t(as.data.frame(result$summary))[1,1:6]
        }
        if(trait_type=="quant"){
                result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="quant", beta=input$beta_gwas,varbeta=input$varbeta, N=quant_trait_num[num_2],snp = input$variant_id),
          dataset2=list(pvalues=input$pval_nominal_sqtl, type="quant", N=celltype_sample[num_3],snp=input$variant_id), MAF=input$maf)
                test_rs[1,4:9]<-t(as.data.frame(result$summary))[1,1:6]}}


if(nrow(test_rs)>0){output_dir <- paste0(celltype, "_", chr, "_nominals_1/", gwas_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(test_rs, file = paste0(output_dir, "/", intron), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)}





