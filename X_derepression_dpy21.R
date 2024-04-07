library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)

cec4_WT <- results(DErun, contrast=c("genotype_treatment","cec4_noauxin","N2_noauxin"))
cec4_WT <- data.frame(cec4_WT)
cec4_WT <- drop_na(cec4_WT)

dpy21_WT <- results(DErun, contrast=c("genotype_treatment","dpy21_noauxin","N2_noauxin"))
dpy21_WT <- data.frame(dpy21_WT)
dpy21_WT <- drop_na(dpy21_WT)

dpy21cec4_WT <- results(DErun, contrast=c("genotype_treatment","dpy21_cec4_noauxin","N2_noauxin"))
dpy21cec4_WT <- data.frame(dpy21cec4_WT)
dpy21cec4_WT <- drop_na(dpy21cec4_WT)

dpy21cec4_cec4 <- results(DErun, contrast=c("genotype_treatment","dpy21_cec4_noauxin","cec4_noauxin"))
dpy21cec4_cec4 <- data.frame(dpy21cec4_cec4)
dpy21cec4_cec4 <- drop_na(dpy21cec4_cec4)

dpy21cec4_dpy21 <- results(DErun, contrast=c("genotype_treatment","dpy21_cec4_noauxin","dpy21_noauxin"))
dpy21cec4_dpy21 <- data.frame(dpy21cec4_dpy21)
dpy21cec4_dpy21 <- drop_na(dpy21cec4_dpy21)

annotate_deseq_res <- function(input_dataframe){
  mart <- useDataset("celegans_gene_ensembl", 
                     useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
  genes <- rownames(input_dataframe)
  gene_list <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
    filters = "ensembl_gene_id",
    values = genes,
    mart = mart, 
    useCache = FALSE)
  gene_list <- data.frame(gene_list)
  out_dataframe <- merge.data.frame(input_dataframe, gene_list, by.x= 0 , by.y = "ensembl_gene_id")
  colnames(out_dataframe)[1] <- "EnsemblID"
  return(out_dataframe)
}

unannotated_DE <- list("cec4_WT" = cec4_WT,
                       "dpy21_WT" = dpy21_WT,
                       "dpy21cec4_WT" = dpy21cec4_WT,
                       "dpy21cec4_cec4" = dpy21cec4_cec4,
                       "dpy21cec4_dpy21" = dpy21cec4_dpy21)

annotated_DE <- lapply(unannotated_DE, annotate_deseq_res)

annotate_XorA <- function(input_dataframe){
  XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
  for (i in 1:nrow(input_dataframe)) {
    if (input_dataframe[i,"chromosome_name"] == "X"){
      XorA[i,1] <- "X chromosome"
    } else {
      XorA[i,1] <- "Autosome"
    }
  }
  XorA$EnsemblID <- input_dataframe$EnsemblID
  out_dataframe <- merge.data.frame(input_dataframe, XorA, by = "EnsemblID")
  return(out_dataframe)
}

annotated_DE_XA <- lapply(annotated_DE, annotate_XorA) 
annotated_DE_XA <- lapply(annotated_DE_XA, subset, chromosome_name != "MtDNA")

assign_factor_level <- function(inputdataframe) {
  inputdataframe$XorA <- as.factor(inputdataframe$XorA)
  inputdataframe$XorA <- relevel(inputdataframe$XorA, "X chromosome")
  return(inputdataframe)
}

annotated_DE_XA <- lapply(annotated_DE_XA, assign_factor_level)

#boxplots
dir.create("boxplots")

ggplot(subset(annotated_DE_XA[["cec4_WT"]], baseMean > 1), 
       aes(x = XorA, y = log2FoldChange, fill = XorA)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#007991","lightgrey")) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  geom_hline(yintercept = median(annotated_DE_XA[["cec4_WT"]]$log2FoldChange[annotated_DE_XA[["cec4_WT"]]$XorA == "Autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 20), 
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(y = "log2 Fold Change", title = "cec-4 vs WT")

ggsave("boxplots/cec4_WT.png", 
       height = 4, width = 5)

ggplot(subset(annotated_DE_XA[["dpy21_WT"]], baseMean > 1), 
       aes(x = XorA, y = log2FoldChange, fill = XorA)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#439A86","lightgrey")) +
  coord_cartesian(ylim = c(-2.5,3.5)) +
  geom_hline(yintercept = median(annotated_DE_XA[["dpy21_WT"]]$log2FoldChange[annotated_DE_XA[["dpy21_WT"]]$XorA == "Autosome"]), 
             linetype = 2, alpha = 0.6) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 20), 
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(y = "log2 Fold Change", title = "dpy-21 vs WT")

ggsave("boxplots/dpy21_WT.png", 
       height = 4, width = 5)

ggplot(subset(annotated_DE_XA[["dpy21cec4_WT"]], baseMean > 1), 
       aes(x = XorA, y = log2FoldChange, fill = XorA)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#BCD8C1","lightgrey")) +
  coord_cartesian(ylim = c(-2.3,3.3)) +
  geom_hline(yintercept = median(annotated_DE_XA[["dpy21cec4_WT"]]$log2FoldChange[annotated_DE_XA[["dpy21cec4_WT"]]$XorA == "Autosome"]), 
             linetype = 2, alpha = 0.6) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 20), 
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(y = "log2 Fold Change", title = "dpy-21; cec-4 vs WT")

ggsave("boxplots/dpy21cec4_WT.png", 
       height = 4, width = 5)

ggplot(subset(annotated_DE_XA[["dpy21cec4_cec4"]], baseMean > 1), aes(x = XorA, y = log2FoldChange, fill = XorA)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#BCD8C1","lightgrey")) +
  coord_cartesian(ylim = c(-2.1,2.7)) +
  geom_hline(yintercept = median(annotated_DE_XA[["dpy21cec4_cec4"]]$log2FoldChange[annotated_DE_XA[["dpy21cec4_cec4"]]$XorA == "Autosome"]), 
             linetype = 2, alpha = 0.6) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 20), 
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(y = "log2 Fold Change", title = "dpy-21; cec-4 vs cec-4")

ggsave("boxplots/dpy21cec4_cec4.png", 
       height = 4, width = 5)

ggplot(subset(annotated_DE_XA[["dpy21cec4_dpy21"]], baseMean > 1), 
       aes(x = XorA, y = log2FoldChange, fill = XorA)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#BCD8C1","lightgrey")) +
  coord_cartesian(ylim = c(-1.6,1.5)) +
  geom_hline(yintercept = median(annotated_DE_XA[["dpy21cec4_dpy21"]]$log2FoldChange[annotated_DE_XA[["dpy21cec4_dpy21"]]$XorA == "Autosome"]), 
             linetype = 2, alpha = 0.6) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 20), 
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(y = "log2 Fold Change", title = "dpy-21; cec-4 vs dpy-21")

ggsave("boxplots/dpy21cec4_dpy21.png", 
       height = 4, width = 5)

##statistical test

run_wilcox_test <- function(inputdataframe){
  p_val <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$XorA == "Autosome" & inputdataframe$baseMean > 1],
                      inputdataframe$log2FoldChange[inputdataframe$XorA == "X chromosome" & inputdataframe$baseMean > 1])[3]
  return(p_val)
  }

stat_XA <- lapply(annotated_DE_XA, run_wilcox_test)
