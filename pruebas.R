path_counts <- "/home/pablo/Documents/practicas/aln/counts.txt"
counts <- read.csv(path_counts, sep="\t", skip = 1)
counts <- counts[,c(1,7:13)]
counts_df <- as.data.frame(counts)

new_colnames <- c("untreated_1", "untreated_2", "untreated_3", "dexa_1", "dexa_2", "dexa_3", "dexa_4")
colnames(counts_df)[2:8] <- new_colnames

rownames(counts_df) <- counts_df[,1]
counts_df <- counts_df[,-1]

annotations_df <- data.frame(
  condition = c("untreated", "untreated", "untreated", "treated", "treated", "treated", "treated"),
  type = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))

rownames(annotations_df) <- c("untreated_1", "untreated_2", "untreated_3", "dexa_1", "dexa_2", "dexa_3", "dexa_4")
annotations_df$condition = factor(annotations_df$condition)
annotations_df$type = factor(annotations_df$type)

library("DESeq2")
deseq_ds <- DESeqDataSetFromMatrix(countData = counts_df, colData = annotations_df, design = ~ condition)
smallestGroupSize <- 3 # 3 untreated y 4 treated
keep <- rowSums(counts(deseq_ds) >= 10) >= smallestGroupSize
deseq_ds <- deseq_ds[keep,]
deseq_ds <- DESeq(deseq_ds)
results <- results(deseq_ds)
coef <- resultsNames(deseq_ds)[2] # coef == "condition_untreated_vs_treated"
library("apeglm")
results_LFC <- lfcShrink(deseq_ds, coef=coef, type="apeglm")

results_df <- as.data.frame(results_LFC)
results_df$gene_id <- rownames(results_df)

library(EnsDb.Hsapiens.v86)
ensemblsIDS <- rownames(results_df)
length(ensemblsIDS) == length(unique((ensemblsIDS)))
columns(EnsDb.Hsapiens.v86)
mapIds <- ensembldb::select(EnsDb.Hsapiens.v86, keys = ensemblsIDS, keytype = "GENEID", column = c("SYMBOL", "GENEID", "TXBIOTYPE"))
mapIds <- mapIds[mapIds$TXBIOTYPE == "protein_coding",]
duplicated_symbols <- mapIds[duplicated(mapIds$SYMBOL) | duplicated(mapIds$SYMBOL, fromLast = TRUE), ]
duplicated_ensg <- duplicated_symbols[,2]
results_df <- results_df[!results_df$gene_id %in% duplicated_ensg,]
rownames(results_df) <- ifelse(results_df$gene_id %in% mapIds$GENEID, mapIds$SYMBOL[match(results_df$gene_id, mapIds$GENEID)], results_df$gene_id)

library(EnhancedVolcano)
EnhancedVolcano(results_df,
                lab = rownames(results_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.00001,
                FCcutoff = 2)



