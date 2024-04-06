library(tidyr)
library(AnnotationHub)
library(ensembldb)
library(tximport)
library(DESeq2)
library(stringr)

deseq2_metadata <- data.frame(sampleID = rep("PH", 33))
sampleIDs <- list.files("salmon_counts/")
sampleIDs <- str_replace(sampleIDs, "salmon_counts_", "")
deseq2_metadata$sampleID <- sampleIDs

metadata_raw <- read.table("Trombley_2024/metadata.csv", 
                           header = TRUE, 
                           sep = ",")
deseq2_metadata <- merge.data.frame(deseq2_metadata, 
                                    metadata_raw[,c("Sample.Name", "File.Designation")], 
                                    by.x = "sampleID", by.y = "File.Designation")
deseq2_metadata$genotype <- deseq2_metadata$Sample.Name
deseq2_metadata$genotype <- deseq2_metadata$genotype %>% 
  str_replace("_rep[123]","") %>% 
  str_replace("_aux", "") %>% 
  str_replace("_noaux", "")

deseq2_metadata$auxin <- rep("noauxin", nrow(deseq2_metadata))
for (i in 1:nrow(deseq2_metadata)) {
  if (deseq2_metadata$Sample.Name[i] %like% "_aux") {
    deseq2_metadata$auxin[i] <- "auxin"
  }
}

deseq2_metadata$genotype_treatment <- paste0(deseq2_metadata$genotype, "_", deseq2_metadata$auxin)
deseq2_metadata$genotype_treatment <- as.factor(deseq2_metadata$genotype_treatment)
summary(deseq2_metadata$genotype_treatment)
colnames(deseq2_metadata)[2] <- "genotype_treatment_rep"

directory <- "salmon_counts/"
list.files(directory)
files <- file.path(directory, 
                   paste0("salmon_counts_", deseq2_metadata$sampleID), 
                   "quant.sf")
names(files) <- deseq2_metadata$sampleID
files
all(file.exists(files))

celegansdb_formalclassobject <- query(AnnotationHub(), 
                                      pattern = c("Caenorhabditis elegans", "EnsDb", "109"))
celegansdb <- celegansdb_formalclassobject[["AH109526"]]
genes <- genes(celegansdb)
tx2gene <- data.frame(TXNAME=genes$canonical_transcript,
                      GENEID=genes$gene_id)
head(tx2gene)

salmon_import <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_import)
head(salmon_import$counts)

dds <- DESeqDataSetFromTximport(salmon_import,
                                colData = deseq2_metadata,
                                design = ~ genotype_treatment)
dds$genotype_treatment <- relevel(dds$genotype_treatment, ref = "N2_noauxin")
DErun <- DESeq(dds)
DErun_res <- results(DErun)
head(DErun_res)

#Save RDS
saveRDS(DErun,"Trombley_2024/DE_run.rds")

rld2 <- vst(DErun, blind = TRUE)
plotPCA(rld2, intgroup = "genotype_treatment", ntop = 100)