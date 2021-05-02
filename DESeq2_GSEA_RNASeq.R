#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
if (!require(DESeq2)) BiocManager::install("DESeq2")


library(TCGAbiolinks) # load TCGAbiolinks library
library(SummarizedExperiment)  # load SummarizedExperiment library
library(DESeq2) # load DESeq2 library

# load GDC data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code once to download the data
sum_exp <- GDCprepare(query)  # create sum_exp variable



### RNA Expression ###

patient_data <- colData(sum_exp)  # save colData in patient_data variable
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis  # set patient_ages as age of initial diagnosis
patient_data$age_category = ifelse(patient_ages < 50, "Young", "Old")  # create new age_category column
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"  # get gene expression counts

patient_data$ESR1_counts = htseq_counts["ENSG00000091831",]  # create column of ESR1 counts in patient_data

patient_data$ESR1_counts_log = sapply(htseq_counts["ENSG00000091831",], log10) # create column of log(ESR1) counts in patient_data

all_young <- sum(patient_data$age_category == "Young", na.rm = TRUE)  # get # of young patients
all_old <- sum(patient_data$age_category == "Old", na.rm=TRUE)  # get # of old patients

patient_data <- subset(patient_data, patient_data$paper_BRCA_Subtype_PAM50 == "LumB" | patient_data$paper_BRCA_Subtype_PAM50 == "LumA")  # create subset of data with only ER-positive
ER_young <- sum(patient_data$age_category == "Young", na.rm = TRUE)  # get # of er-positive young patients
ER_old <- sum(patient_data$age_category == "Old", na.rm=TRUE) # get # of er-positive old patients

ER_young/all_young  # get percentage of young patients that are ER-positive
ER_old/all_old  # get percentage of old patients that are ER-positive

require(stats)  # load stats
jpeg("log(ESR1)_RNA_scatterplot")  # create jpeg
reg_log <- lm(ESR1_counts_log ~ paper_age_at_initial_pathologic_diagnosis, data = patient_data)  # perform regression on log(ESR1)
coeff_log <- coefficients(reg_log)  # get coefficients to get equation
coeff_log  # print coefficients
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$ESR1_counts_log, main = "y = 0.01x + 3.87")  # plot scatterplot
abline(reg_log)  # plot regression line
dev.off()  # plot complete


jpeg("ESR1_RNA_scatterplot")  # create jpeg
reg <- lm(ESR1_counts ~ paper_age_at_initial_pathologic_diagnosis, data = patient_data)  # regression on ESR1
coeff <- coefficients(reg)  # get coefficients to get equation
coeff  # print coefficients
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$ESR1_counts, main = "y = 1015.7x - 14201.4")  # plot scatterplot
abline(reg)  # plot regression line 
dev.off()  # plot complete


jpeg("log(ESR1)_RNA_boxplot")  # create jpeg
boxplot(ESR1_counts_log~factor(age_category, levels=c("Young", "Old")), data = patient_data, main = "Boxplot of log(HTSeq - Counts) for ESR1 by Age Category")  # plot boxplot of log(ESR1)
dev.off()  # plot complete

library(ggpubr)  # load ggpubr library
patient_data$age_category_factor <- as.factor(patient_data$age_category)  # make age_category a factor
patient_data <- as.data.frame(patient_data) # make sure patient_data is a dataframe
compare_means(ESR1_counts ~ age_category, data=patient_data)  # calculates p-value



### DESeq2 ###

#access counts
counts <- assays(sum_exp)$"HTSeq - Counts"
# remove patients with NA in age, PAM50, pathology
patients_no_NA_mask <- ( !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis)
                         & !is.na(colData(sum_exp)$paper_BRCA_Subtype_PAM50)
                         & !is.na(colData(sum_exp)$paper_BRCA_Pathology)
                         & !colData(sum_exp)$paper_BRCA_Pathology == "NA" )
#access the patient_data from coldata and apply NA mask
patient_data <- colData(sum_exp)[ patients_no_NA_mask, ]

# apply NA mask to counts
counts <- counts[ , patients_no_NA_mask]

# create age_category column
patient_data$age_category <- ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 50, "Young", "Old")


# make categorical columns factors 
patient_data$age_category <- factor(patient_data$age_category, levels=c("Young", "Old"))
patient_data$paper_BRCA_Subtype_PAM50 <- factor( patient_data$paper_BRCA_Subtype_PAM50, levels=c("Her2","LumA","LumB","Basal","Normal") )
patient_data$paper_BRCA_Pathology <- factor( patient_data$paper_BRCA_Pathology, levels=c("IDC","Other","Mixed","ILC") )

# prepare DESeq2 data matrix
dds_with_adjustment <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~paper_BRCA_Pathology+ paper_BRCA_Subtype_PAM50 +age_category)
# perform DESeq2 analysis
dds_obj <- DESeq(dds_with_adjustment)
# load results 
results <- results(dds_obj, contrast=c("age_category", "Young",'Old'))
# create mask of NA values
results_NA <- !is.na(results[,"padj"])
# apply mask to results
results <- results[results_NA,]
# create new column of fold change
results$FoldChange <- 2^results$log2FoldChange
# save results csv file
write.csv(results, "/PATHWAY/DESeq2.csv")
# separate significant genes from insignificant
results_significant_adjp <- results[results$padj < padj_threshold,]  #separate significant genes from insignificant ones
# separate up regulated in young
results_sig_up_regulated <- results_significant_adjp[results_significant_adjp$log2FoldChange > log2FC_threshold,]  
# separate down regulated in young
results_sig_down_regulated <- results_significant_adjp[results_significant_adjp$log2FoldChange < -log2FC_threshold,]

gene_information <- rowData(sum_exp)  #access gene information from sum_exp rowData
matching_rows_up <- gene_information[gene_information$ensembl_gene_id %in% rownames(results_sig_up_regulated),]  # create dataframe with only gene information for genes that were up regulated
results_sig_up_regulated$gene_name <- matching_rows_up$external_gene_name # create new column in up regulated df for common gene name
matching_rows_down <- gene_information[gene_information$ensembl_gene_id %in% rownames(results_sig_down_regulated),] # create dataframe with only gene information for genes that were down regulated
results_sig_down_regulated$gene_name <- matching_rows_down$external_gene_name  # create new column in down regulated df for common gene names

#save as .csv files
write.csv(results_sig_up_regulated, "/PATHWAY/high_in_young.csv")
write.csv(results_sig_down_regulated, "/PATHWAY/low_in_young.csv")


### GSEA .rnk File ###

# create variable of gene names of genes in results
matching_rows_up <- gene_information[gene_information$ensembl_gene_id %in% rownames(results),]  
# create new gene_name column in results 
results$gene_name <- matching_rows_up$external_gene_name 
# create gsea_rnk dataframe with gene names column
gsea_rnk <- subset(results, select=c("gene_name"))
# add log2FC column
gsea_rnk$log2FC <- results$log2FoldChange
# write dataframe to .rnk file
write.table(gsea_rnk, "/PATHWAY/file.rnk", sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)



