

/public/home/liangyan/software/regtools junctions extract -a 5 -m 30 -M 500000 -s XS bam -o bam.junc

python leafcutter_cluster_regtools.py -j test_juncfiles.txt_celltype -m 20 -o oneK1K_celltype -l 500000

python prepare_phenotype_table.py oneK1K_celltype_perind.counts.gz

##intron_filter
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]

dat <- read.delim(paste0('oneK1K_celltype_perind_numers.counts.gz'), sep=' ', check.names=F)
introns <- c()
colNum <- ncol(dat)
for (i in 1:nrow(dat)){
  if (length(which(dat[i,] == 0))/colNum <= 0.6){
    introns <- c(introns, rownames(dat)[i])
  }
}
write.table(introns, paste0('filtered_celltype.txt'), col.names=F, row.names=F, quote=F)


##Demographic-biased AS
Rscript_path="/public/home/liangyan/software/R/bin/Rscript"
script_path="/public/home/liangyan/software/leafcutter/scripts/leafcutter_ds.R"
output_dir="./result/"
input_dir="./junc/"
ct_dir="./ct/"

for celltype in "${celltypes[@]}"; do
    echo "Processing $celltype..."
    output="${output_dir}leafcutter_ds.sex.${celltype}"
    junc_file="${input_dir}S20_${celltype}/oneK1K_${celltype}_perind_numers.counts.gz"
    ct_file="${ct_dir}${celltype}.txt"

    $Rscript_path $script_path -o $output $junc_file $ct_file
done



