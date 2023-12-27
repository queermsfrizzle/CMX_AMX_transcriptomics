install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")

library(dada2)

pathF <- "~/scratch/AJP10_Amplicon_Seq/16S/forward/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "~/scratch/AJP10_Amplicon_Seq/16S/reverse/" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
trimoutput <-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,240), truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

head(trimoutput)


filtpathF <- "~/scratch/AJP10_Amplicon_Seq/16S/forward/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "~/scratch/AJP10_Amplicon_Seq/16S/reverse/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "~/scratch/AJP10_Amplicon_Seq/16S/seqtabmidas.rds") # CHANGE ME to where you want sequence table saved


st.all <- readRDS("~/scratch/AJP10_Amplicon_Seq/16S/seqtabmidas.rds")
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/storage/scratch1/9/jjohnston74/AJP10_Amplicon_Seq/16S/DADA2 File MiDAS 5.1.fa", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "~/scratch/AJP10_Amplicon_Seq/16S/seqtab_midas.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "~/scratch/AJP10_Amplicon_Seq/16S/tax_final_midas.rds") # CHANGE ME ...

write.csv(seqtab, "~/scratch/AJP10_Amplicon_Seq/16S/seqtab_final_midas.csv", col.names=NA)
write.csv(tax, "~/scratch/AJP10_Amplicon_Seq/16S/tax_final_midas.csv", col.names=NA)

seqtab.flip <-t(seqtab)
write.csv(seqtab.flip, "~/scratch/AJP10_Amplicon_Seq/16S/seqtab_transposed_final_midas.csv", col.names=NA)

---------------------------------
  
library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)

nproc <- 4 
asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters",
  processors = nproc)

clusters <- clusters %>%
  add_column(sequence = asv_sequences)
merged_seqtab <- seqtab %>%
  # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
  t %>%
  as_tibble(rownames = "sequence") %>%
  # add the cluster information
  left_join(clusters, by = "sequence") %>%
  # merge ASVs in the same cluster, summing abundances within samples
  group_by(cluster) %>%
  summarize_at(vars(-sequence), sum) %>%
  # Set new taxa names to OTU<cluster #> 
  mutate(cluster = paste0("OTU", cluster)) %>%
  # convert back to a matrix in the original orientation
  column_to_rownames("cluster") %>%
  as("matrix") %>%
  t

otu_16s_seqtab <-t(merged_seqtab)
write.csv(otu_16s_seqtab, "~/scratch/AJP10_Amplicon_Seq/16S/otu_16s_seqtab_midas.csv", col.names=NA)
write.csv(clusters, "~/scratch/AJP10_Amplicon_Seq/16S/clusters_midas.csv", col.names=NA)

