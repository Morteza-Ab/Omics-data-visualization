setwd("/Users/mortezaabyadeh/Desktop")
library(data.table)
library(ggplot2)
library(DESeq2)
library(openxlsx) #don't need here but in case for other analysis I kept it here
library(readxl)
data <- read_excel("combined_proteome_transcriptome_data.xlsx")
class(data)
colnames(data)

data[,1] <- rownames(data)
data <- data[,-1]
gr <- c(rep("proteome-lsm5", 3), rep("proteome-ctrl", 3), rep("transcriptome-lsm5", 3), rep("transcriptome-ctrl", 3))
gr <- as.factor(gr)

########we had float numbers that changes ti integer

data <- data[, sapply(data, is.numeric)] <- lapply(data[, sapply(data, is.numeric)], as.integer)

### I had NA in the data: if (any(is.na(data))) {
# print("There are NA values in the count matrix.")

## total_na <- sum(is.na(data))
## print(paste("Total number of NA values:", total_na)) 

cleaned_data <- na.omit(data)

colData <- data.frame(Group = gr)
rownames(colData) <- colnames(cleaned_data)
cds <- DESeqDataSetFromMatrix(countData = cleaned_data, colData = colData, design = ~ Group)
cds <- DESeq(cds)
head(data)


data.norm <- log2(1+counts(cds, normalized=T))
head(data.norm)
boxplot(data.norm)

###Function_I made that_can help to more discriminate between data in pca plot-You may do not like it!
data.mean.center <- t(scale(t(data.norm), scale = F))


###pcr
pc <- prcomp(data.mean.center)
plot(pc)
head(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
head(data)
pcr$group <- gr
head(pcr)


library(ggplot2)

ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 5, alpha = 0.9) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=5, alpha=0.9, shape = 10) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                                                         axis.text = element_text(size = 12),
                                                                                                         axis.title = element_text(size = 14, face = "bold"),
                                                                                                         legend.title = element_text(size = 14, face = "bold"),
                                                                                                         legend.text = element_text(size = 12),
                                                                                                         legend.position = "right")
cbPalette <- c("#999999", "#E69F00", "#56B4E9")

custom_palette <- c("lightgreen", "gold", "lightblue","purple")

ggplot(pcr, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 9, alpha = 0.9, position = position_jitter(width = 0.1, height = 0.1)) +  
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "right")
