BiocManager::install(c("scater","scran"))
library(scran)
library(scater)
library(Matrix)
synapser::synLogin()

p <- synapser::synGet('syn18686381')
counts <- readMM(p$path)

# get sample QC/batch data: snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv
syn18642934 <- synapser::synGet(entity='syn18642934') 
batches <- read.csv(syn18642934$path,stringsAsFactors = F)

#get tags and celltype data: filtered_column_metadata.txt
p2 <- synapser::synGet('syn18686372')
Labels <- read.delim(p2$path,stringsAsFactors = F)

#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')
rownames(counts) <- readLines(p3$path)   



# get rosmap id mappings: snRNAseqPFC_BA10_id_mapping.csv
p3 <- synapser::synGet('syn18694015')
ids <- read.csv(p3$path,stringsAsFactors = F)

####upload a csv version of supplemtary table 5 from mathys et al. into R to get neuropath data
metadata <- read.csv('~/scRNAseq-workshop2020/Mathys_supplement_5.csv')
#ROSMAP_clinical.csv
p4 <- synapser::synGet(entity='syn3191087') 
metadata2 <- read.csv(p4$path)

#extract the mathys patients only and impute the missing PMI as the median PMI of the dataset, add that to 
#patient missing pmi in metadata file. Also impute missing apoe genotypes to 3/3
metadata3 <- metadata2
metadata3 <- plyr::match_df(metadata3, ids, on="projid")

paste('Imputing PMI to:',median(metadata3$pmi[!is.na(metadata3$pmi)]))
#add this back into the metadata file
metadata2$pmi[is.na(metadata2$pmi)] <- 7
#impute missing apoe genotypes as '33'
metadata2$apoe_genotype[is.na(metadata2$apoe_genotype)] = 33
head(metadata2)


#create list categories for the features data frame
sex = c()
m = as.character(metadata$msex)
fileName = c()
batch = c()
ros_ids = c()
projid = c()
Diagnosis = c()
cogdx = c()
ceradsc = c()
braaksc = c()
tangles = c()
apoe_genotype = c()
pmi = c()
educ=c()
race=c()
AOD=c()

#Need to be able to harmonize labels/counts/metadata by the different identifiers (projid, Subject, rosids)
for(i in 1:length(rownames(Labels))){
  ros_ids = c(ros_ids,ids$Subject[c(which(ids$projid==Labels$projid[i])[1])])
  fileName = c(fileName,ids$fastq[c(which(ids$projid==Labels$projid[i])[1])])
}
Labels$ros_ids = ros_ids
Labels$fileName = fileName
head(Labels)

#map metadata to cells
for(i in 1:length(rownames(Labels))){
  batch = c(batch,batches$sequencingBatch[c(which(batches$fileName==Labels$fileName[i])[1])])
  Diagnosis = c(Diagnosis,metadata$pathology.group[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  cogdx = c(cogdx,metadata$cogdx[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  sex = c(sex,m[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  ceradsc = c(ceradsc,metadata$ceradsc[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  braaksc = c(braaksc,metadata$braaksc[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  tangles = c(tangles,metadata$tangles[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  apoe_genotype = c(apoe_genotype,metadata2$apoe_genotype[c(which(metadata2$projid==Labels$projid[i])[1])])
  pmi = c(pmi,metadata2$pmi[c(which(metadata2$projid==Labels$projid[i])[1])])
  educ = c(educ,metadata2$educ[c(which(metadata2$projid==Labels$projid[i])[1])])
  race = c(race,metadata2$race[c(which(metadata2$projid==Labels$projid[i])[1])])
  AOD = c(AOD,metadata2$age_first_ad_dx[c(which(metadata2$projid==Labels$projid[i])[1])])
}


Labels$batch = batch
Labels$Diagnosis = Diagnosis
Labels$cogdx = cogdx
Labels$sex = sex
Labels$ceradsc = ceradsc
Labels$braaksc = braaksc
Labels$tangles = tangles
Labels$apoe_genotype = apoe_genotype
Labels$pmi = pmi
Labels$educ=educ
Labels$race=race
Labels$AOD=AOD

Labels$Diagnosis[Labels$Diagnosis=='late-pathology'] <- "Late-pathology AD"
Labels$Diagnosis[Labels$Diagnosis=='no-pathology'] <- 'Control'
Labels$Diagnosis[Labels$Diagnosis=='early-pathology'] <- 'Early-pathology AD'
Labels$simpleDiagnosis = Labels$Diagnosis
Labels$simpleDiagnosis[Labels$simpleDiagnosis!='Control'] <- "AD"
Labels$Diagnosis[Labels$Diagnosis==1]='Early-pathology AD'
Labels$Diagnosis[Labels$Diagnosis==2]='Late-pathology AD'
Labels$Diagnosis[Labels$Diagnosis==3]='Control'
head(Labels)

colnames(counts) <- Labels[,1]
rownames(Labels) <- Labels[,1]

gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
head(gene_short_name)
length(gene_short_name)



sce <- SingleCellExperiment(list(counts=counts))
dim(sce)
clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)

summary(sizeFactors(sce))

sce <- logNormCounts(sce)
head(logcounts(sce[0:20,0:20]))
head(rownames(logcounts(sce)))

dim(logcounts(sce))

counts2 = logcounts(sce)

saveRDS(sce, file="mathys_scran_normalized.rds")
saveRDS(counts2, file="mathys_scran_normalized_counts2.rds")



