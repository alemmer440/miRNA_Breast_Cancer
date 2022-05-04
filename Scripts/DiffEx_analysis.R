#Differential Expression Analysis of TCGA-BRCA Data- miRNAs

library( "DESeq2" )
library(ggplot2)

#Read in miRNA counts data and phenotype files
miRNA_counts<- read.csv('../Organized_Data/miRNA_counts_v2.csv', row.names=1, header=TRUE)

miRNA_pheno<- read.csv('../Organized_Data/miRNA_pheno_v2.csv', header=TRUE)

#Use DESeq to conduct diff expression analysis 
dds <- DESeqDataSetFromMatrix(countData = miRNA_counts,
                              colData = miRNA_pheno,
                              design= ~Phenotype)

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)

head(results(dds, tidy=TRUE))

#order by p values
res <- res[order(res$padj),]
head(res)

#Save results
write.csv(as.data.frame(res), 
         file="miRNA_DESea_results.csv")


#Label miRNAs as signficant or not and upexpressed or downexpressed
resultdf=as.data.frame(res)
# add a column of NAs
resultdf$diffexpressed <- "Not Significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
resultdf$diffexpressed[resultdf$log2FoldChange > 1.5 & resultdf$padj < 0.05] <- "Significant"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
resultdf$diffexpressed[resultdf$log2FoldChange < -1.5 & resultdf$padj < 0.05] <- "Significant"
#write.csv(as.data.frame(resultdf), 
#  file="miRNA_DESeq_results_dfv2.csv")

p <- ggplot(data=resultdf, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + geom_point(size=4) +  scale_color_manual(values = c("Not Significant" = "black",
                                                                                                                                     "Significant"="red")) +
  ggtitle("Differential miRNA expression in breast cancer \ntumors compared with healthy tissue")+ theme(plot.title = element_text(hjust = 0.5))+ ylab("-log10(pvalue)")+ theme(legend.title = element_blank())+ theme(text = element_text(size = 35))
#p<-p+ ggtitle("Tumor vs Normal Sample")+ ylab("-log10(pvalue")
#p<-p + theme(legend.title = element_blank())
ggsave('volcanoplotv12.jpg', plot=p, width = 15, height = 15, dpi = 150, units = "in")

#sig_miRNA= resultdf[resultdf$diffexpressed=='Significant',]
#write.csv(as.data.frame(sig_miRNA), file='significant_miRNAs.csv')
          
