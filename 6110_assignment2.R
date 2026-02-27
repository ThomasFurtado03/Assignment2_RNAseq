#Assignment 2
#Thomas Furtado
sessionInfo()

#Necessary Packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
install.packages("pheatmap")
install.packages("ggplot2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")

#Load Packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Sc.sgd.db)




samples <- data.frame(
  sample = c("SRR10551665","SRR10551664","SRR10551663",
             "SRR10551662","SRR10551661","SRR10551660",
             "SRR10551659","SRR10551658","SRR10551657"),
  condition = factor(c("Early", "Early", "Early",
                       "Thin", "Thin", "Thin",
                       "Mature", "Mature", "Mature"),
                     levels = c("Early", "Thin", "Mature"))
)

rownames(samples) <- samples$sample
samples

files <- file.path("kallisto_output",
                   samples$sample,
                   "abundance.tsv")
names(files) <- samples$sample

files

txi <- tximport(files, type = "kallisto", txOut = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

dds
resultsNames(dds)

#Contrasts
res_Thin_vs_Early   <- results(dds, name = "condition_Thin_vs_Early")
res_Mature_vs_Early <- results(dds, name = "condition_Mature_vs_Early")
res_Mature_vs_Thin  <- results(dds, contrast = c("condition","Mature","Thin"))

summary(res_Thin_vs_Early)
summary(res_Mature_vs_Early)
summary(res_Mature_vs_Thin)


#Export DE tables
write.csv(as.data.frame(res_Thin_vs_Early), "DE_Thin_vs_Early.csv")
write.csv(as.data.frame(res_Mature_vs_Early), "DE_Mature_vs_Early.csv")
write.csv(as.data.frame(res_Mature_vs_Thin),  "DE_Mature_vs_Thin.csv")

exists("res_Thin_vs_Early")
exists("res_Mature_vs_Early")
exists("res_Mature_vs_Thin")


#Overall data structure: PCA

vsd <- vst(dds, blind = FALSE)
PCA <- plotPCA(vsd, intgroup = "condition")
print(PCA)



#Top DE genes

#order by adjusted p-value
res_ordered <- res_Mature_vs_Early[order(res_Mature_vs_Early$padj), ]

#Collect top 20 genes
top_genes <- rownames(res_ordered)[1:20]

#Normalized counts
norm_counts <- assay(vsd)[top_genes, ]

pheatmap(norm_counts, scale="row",
         fontsize_row=6)



#Clean up gene names
gene_names <- rownames(res_Mature_vs_Early)
gene_names <- sub("_mRNA", "", gene_names)
rownames(res_Mature_vs_Early) <- gene_names

head(rownames(res_Mature_vs_Early))


#Define significant genes
res_df <- as.data.frame(res_Mature_vs_Early)

sig_genes <- rownames(res_df)[
  which(!is.na(res_df$padj) & 
          res_df$padj < 0.05 & 
          abs(res_df$log2FoldChange) > 1)
]

length(sig_genes)

#Define gene background
background_genes <- rownames(res_df)


#ORA
keyType = "ORF"

ego <- enrichGO(
  gene = sig_genes,
  universe = background_genes,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
)

nrow(as.data.frame(ego))

write.csv(as.data.frame(ego), "results/GO_ORA_BP_Mature_vs_Early.csv", row.names = FALSE)


#Functional enrichment figure
p <- dotplot(ego, showCategory = 10, title = "GO Biological Processes")
p <- p + theme(axis.text.y = element_text(size = 8))

p

ggsave("results/Figure_GO_Dotplot_Top10.png", plot = p, width = 7, height = 5, dpi = 300)
