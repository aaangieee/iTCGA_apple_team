# # install packages (you only need to do this once)
# install.packages("tidyverse")
# install.packages("ggrepel")
# 
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("edgeR")

# load the packages you need 
library(tidyverse)
library(ggrepel)
library(edgeR)

# read in the data
gene.counts <- read_table("/Users/jiaqiwang/Desktop/lung_cancer_all_counts.txt", skip = 1)

# set the treatments
condition <- c("DMSO", "DMSO", "DMSO", "Mitoxantrone", "Mitoxantrone", "Mitoxantrone")

# filter out genes expressed at a low level across samples
dim(gene.counts)

totalexp <- rowSums(gene.counts[,7:12])

hist(totalexp)

gene.counts <-  filter(gene.counts, totalexp > 10)

# put the annotation information in a separate dataframe
ann <- gene.counts[,1:6]

# modify the data format for edgeR
d <- DGEList(counts=gene.counts[,7:12], group=factor(condition), genes = ann)
str(d)
dim(d) #11,591 genes in 9 samples

# normalize the data
d <- calcNormFactors(d, method ="TMM")
d$samples

# plot MDS
samples <- c("DMSO_S1", "DMSO_S2", "DMSO_S3", "Mitoxantrone_S1", "Mitoxantrone_S2", "Mitoxantrone_S3")

plotMDS(d, col=as.numeric(d$samples$group), labels = d$samples$group)
legend("center", as.character(unique(d$samples$group)), col=1:3, pch=20)

# start to fit a model (this won't work unless you have more replication)
dm <- estimateCommonDisp(d, verbose = T)
dm <- estimateTagwiseDisp(dm)
plotBCV(dm)
# if this were nice and flat (following the red line, we'd be happy with it, but it looks like we need a more complex model)

# generalized linear model fit
design <- model.matrix(~ 0 + d$samples$group)
design
colnames(design) <- levels(d$samples$group)

dg <- estimateGLMCommonDisp(d, design)
dg <- estimateGLMTrendedDisp(dg, design)
dg <- estimateGLMTagwiseDisp(dg, design)
plotBCV(dg)
# this looks like a much better fit

# Let's fit our new model
fit <- glmFit(dg, design)

# this compares group 1 (Cs, 1) to group 2 (Ts, -1), ignoring group 3 (Vs, 0)
# and does a likelihood ratio test for each gene
fitCT <- glmLRT(fit, contrast=c(1, -1))


deCT <- decideTests(fitCT, adjust.method="BH", p.value = 0.001)
deCTtags <- rownames(dg)[as.logical(deCT)]
plotSmear(fitCT, de.tags=deCTtags)
abline(h = c(-2, 2), col = "blue")

# sort out the differentially expressed genes
tabCT <- topTags(fitCT,n=Inf,adjust.method="BH", sort.by = "PValue")$table

# volcano
# make a significance column
tabCT <- tabCT %>% 
  mutate(significance = case_when((FDR < 0.05 & logFC > 1) ~ "Upregulated", (FDR < 0.05 & logFC < -1) ~ "Downregulated", .default = "Not significant"))

# save this file and also one with just the significant genes
write_delim(tabCT, file = "~/Desktop/count/deCT_2.txt", delim = "\t")
write_delim(filter(tabCT, significance != "Not significant"), file = "~/Desktop/count/deCT_2-sig.txt", delim = "\t")

top_genes <- filter(tabCT, -log10(PValue) > 5)

volcano_plot <- ggplot(tabCT, aes(x = logFC, y = -log10(PValue), color = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not significant" = "grey","Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes, aes(label = Geneid), size = 3.5, fontface = 'bold') +
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  theme_classic() +
  theme(legend.position = 'None',
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank()
  ) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value")

volcano_plot

# Save the volcano plot
ggsave("~/Downloads/volcano_CT_sm.png", plot = volcano_plot, width = 7, height = 5, units = "in", device = "png")

# plot DE of specific genes
# for ENSG00000173110
singlecounts <- gene.counts %>%
  filter(Geneid == "ENSG00000173110") %>% 
  select(-c(Geneid, Chr, Start, End, Strand,Length)) %>% 
  t(.)

singleg <- data.frame(cbind(as.character(d$samples$group), singlecounts))
colnames(singleg) <- c("Treatment","Reads")
singleg$Reads <- as.numeric(singleg$Reads)

sg1 <- ggplot(data = singleg, aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HSPA6") +
  theme_bw()

sg1

ggsave("~/Downloads/HSPA6-de.png", plot = sg1, width = 4, height = 3, units = "in", device = "png")

# for ENSG00000225217
singlecounts <- gene.counts %>%
  filter(Geneid == "ENSG00000225217") %>% 
  select(-c(Geneid, Chr, Start, End, Strand,Length)) %>% 
  t(.)

singleg <- data.frame(cbind(as.character(d$samples$group), singlecounts))
colnames(singleg) <- c("Treatment","Reads")
singleg$Reads <- as.numeric(singleg$Reads)

sg2 <- ggplot(data = singleg, aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HSPA7") +
  theme_bw()

sg2
ggsave("~/Downloads/HSPA7_2-de.png", plot = sg2, width = 4, height = 3, units = "in", device = "png")

# for ENSG00000204388
singlecounts <- gene.counts %>%
  filter(Geneid == "ENSG00000204388") %>% 
  select(-c(Geneid, Chr, Start, End, Strand,Length)) %>% 
  t(.)

singleg <- data.frame(cbind(as.character(d$samples$group), singlecounts))
colnames(singleg) <- c("Treatment","Reads")
singleg$Reads <- as.numeric(singleg$Reads)

sg3 <- ggplot(data = singleg, aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HSPA1B") +
  theme_bw()

sg3

ggsave("~/Downloads/HSPA1B2-de.png", plot = sg3, width = 4, height = 3, units = "in", device = "png")

# for ENSG00000100292
singlecounts <- gene.counts %>%
  filter(Geneid == "ENSG00000100292") %>% 
  select(-c(Geneid, Chr, Start, End, Strand,Length)) %>% 
  t(.)

singleg <- data.frame(cbind(as.character(d$samples$group), singlecounts))
colnames(singleg) <- c("Treatment","Reads")
singleg$Reads <- as.numeric(singleg$Reads)

sg4 <- ggplot(data = singleg, aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HMOX1") +
  theme_bw()
sg4
ggsave("~/Downloads/HMOX1-de.png", plot = sg4, width = 4, height = 3, units = "in", device = "png")



