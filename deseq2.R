library(readxl)
library("DESeq2")
library(tidyverse)
library("writexl")
amx_281_deseq <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "amx281_ifas")

amx_281_deseq <-as.matrix(amx_281_deseq)
row.names(amx_281_deseq) <-amx_281_deseq[,1]
amx_281_cts <-as.data.frame(amx_281_deseq[,-1], stringsAsFactors=FALSE)

amx_281_cts <- amx_281_cts %>% mutate_if(is.character, as.numeric)

amx_281_factor <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "sorting")
amx_281_factor <-as.data.frame(amx_281_factor)
row.names(amx_281_factor) <-amx_281_factor[,1]

all(rownames(amx_281_factor) == colnames(amx_281_cts))
dds <- DESeqDataSetFromMatrix(countData = amx_281_cts,
                                    colData = amx_281_factor,
                                    design = ~ DO)
dds

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

res2 <-as.data.frame(res)
res2$annotation <-row.names(res2)
write_xlsx(res2, path="~/Pinto/DO Transcriptomics Project Data/deseq_results/deseq2_amx281.xlsx")










amx_281_deseq <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "amx_281_ss")

amx_281_deseq <-as.matrix(amx_281_deseq)
row.names(amx_281_deseq) <-amx_281_deseq[,1]
amx_281_cts <-as.data.frame(amx_281_deseq[,-1], stringsAsFactors=FALSE)

amx_281_cts <- amx_281_cts %>% mutate_if(is.character, as.numeric)

amx_281_factor <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "sorting")
amx_281_factor <-as.data.frame(amx_281_factor)
row.names(amx_281_factor) <-amx_281_factor[,1]

all(rownames(amx_281_factor) == colnames(amx_281_cts))
dds <- DESeqDataSetFromMatrix(countData = amx_281_cts,
                              colData = amx_281_factor,
                              design = ~ DO)
dds

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 3) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

res2 <-as.data.frame(res)
res2$annotation <-row.names(res2)
write_xlsx(res2, path="~/Pinto/DO Transcriptomics Project Data/deseq_results/deseq2_amx281_ss_v3.xlsx")





amx_281_deseq <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "cmx_1")

amx_281_deseq <-as.matrix(amx_281_deseq)
row.names(amx_281_deseq) <-amx_281_deseq[,1]
amx_281_cts <-as.data.frame(amx_281_deseq[,-1], stringsAsFactors=FALSE)

amx_281_cts <- amx_281_cts %>% mutate_if(is.character, as.numeric)

amx_281_factor <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "sorting")
amx_281_factor <-as.data.frame(amx_281_factor)
row.names(amx_281_factor) <-amx_281_factor[,1]

all(rownames(amx_281_factor) == colnames(amx_281_cts))
dds <- DESeqDataSetFromMatrix(countData = amx_281_cts,
                              colData = amx_281_factor,
                              design = ~ DO)
dds

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

res2 <-as.data.frame(res)
res2$annotation <-row.names(res2)
write_xlsx(res2, path="~/Pinto/DO Transcriptomics Project Data/deseq_results/deseq2_cmx_1_v3.xlsx")



amx_281_deseq <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "cmx_2")

amx_281_deseq <-as.matrix(amx_281_deseq)
row.names(amx_281_deseq) <-amx_281_deseq[,1]
amx_281_cts <-as.data.frame(amx_281_deseq[,-1], stringsAsFactors=FALSE)

amx_281_cts <- amx_281_cts %>% mutate_if(is.character, as.numeric)

amx_281_factor <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "sorting")
amx_281_factor <-as.data.frame(amx_281_factor)
row.names(amx_281_factor) <-amx_281_factor[,1]

all(rownames(amx_281_factor) == colnames(amx_281_cts))
dds <- DESeqDataSetFromMatrix(countData = amx_281_cts,
                              colData = amx_281_factor,
                              design = ~ DO)
dds

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

res2 <-as.data.frame(res)
res2$annotation <-row.names(res2)
write_xlsx(res2, path="~/Pinto/DO Transcriptomics Project Data/deseq_results/deseq2_cmx_2_v3.xlsx")