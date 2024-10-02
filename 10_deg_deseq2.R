library(DESeq2)
library(rhdf5)
library(tximport)
library(tidyverse)
library(biomaRt)
library(tidyr)
library(dplyr)
# Establishing a biomaRt object/data using the homo sapiens gene ensembl database version 104
#biomart is weird, use this first
httr::set_config(httr::config(ssl_verifypeer = FALSE))
#it used to work without host before
mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=104,
                                                      host = "https://useast.ensembl.org"))
# Creating a list of gene attributes that can be collected in this gene list
attributes_in_mart <- listAttributes(mart)
# build a gene list of all of the transcipt_id, transcript_id_version, ensembl_gene_id, and the hgnc_id values (gene symbols)
# get tss, chr, start and stop site info
# get the new (to v104) canonical transcript information
ensembl_gene_list <- getBM(attributes = c("ensembl_transcript_id","ensembl_transcript_id_version",
											"ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
											"end_position","strand","transcript_start","transcript_end","transcription_start_site",
											"transcript_is_canonical"),mart=mart)
# build the file structure needed by DESeq2/tximport. This is transcriptID then geneID
# select the ensembl_transcript_id and ensembl_gene_id cols and subset it to a new file using
transcript2gene <- ensembl_gene_list %>% dplyr::select(ensembl_transcript_id,ensembl_gene_id)
#get sample annotations and paths to abundance files from a table
csv_mcfkd  = read.csv('file_order_mcf_kd.csv', header = FALSE)
csv_scckd  = read.csv('file_order_scc_kd.csv', header = FALSE)
csv_koscc  = read.csv('file_order_mcfko_scc.csv', header = FALSE)
files_mcfkd <- file.path(csv_mcfkd$V2)
files_scckd <- file.path(csv_scckd$V2)
files_koscc <- file.path(csv_koscc$V2)
names(files_mcfkd) <- paste0(csv_mcfkd$V1)
names(files_scckd) <- paste0(csv_scckd$V1)
names(files_koscc) <- paste0(csv_koscc$V1)
#import kallisto abundance files
txi.mcfkd.kallisto <- tximport(files_mcfkd,type="kallisto", tx2gene=transcript2gene,
                             txOut=FALSE, ignoreTxVersion = TRUE)
txi.scckd.kallisto <- tximport(files_scckd,type="kallisto", tx2gene=transcript2gene,
                             txOut=FALSE, ignoreTxVersion = TRUE)
txi.koscc.kallisto <- tximport(files_koscc,type="kallisto", tx2gene=transcript2gene,
                             txOut=FALSE, ignoreTxVersion = TRUE)
# build a DESeq2 sample table
sampleTable_mcfkd <- data.frame(genotype=factor(c('MCF10A_p63sh','MCF10A_p63sh',
                                                  'MCF10A_SCR','MCF10A_SCR')))
sampleTable_scckd<- data.frame(genotype=factor(c("SCC25_ctr_shRNA","SCC25_ctr_shRNA",
                                                 "SCC25_p63_shRNA","SCC25_p63_shRNA")))
sampleTable_koscc <- data.frame(genotype=factor(c("MCF10Ap53KO_DMSO","MCF10Ap53KO_DMSO","MCF10Ap53KO_DMSO",
                                                  "SCC25_ctr_shRNA","SCC25_ctr_shRNA")))
# row names to each of the value  created in the above data.frame function.
# structure (row name, header and then sample type/condition) minimum for DEseq2
rownames(sampleTable_mcfkd) <- colnames(txi.mcfkd.kallisto$counts)
rownames(sampleTable_scckd) <- colnames(txi.scckd.kallisto$counts)
rownames(sampleTable_koscc) <- colnames(txi.koscc.kallisto$counts)
# creating the DEseq2 object
dds_mcfkd <- DESeqDataSetFromTximport(txi.mcfkd.kallisto, sampleTable_mcfkd, ~genotype)
dds_scckd <- DESeqDataSetFromTximport(txi.scckd.kallisto, sampleTable_scckd, ~genotype)
dds_koscc <- DESeqDataSetFromTximport(txi.koscc.kallisto, sampleTable_koscc, ~genotype)
# performing DESeq2 differential gene expression
mcfkd.de <- DESeq(dds_mcfkd)
scckd.de <- DESeq(dds_scckd)
koscc.de <- DESeq(dds_koscc)
# Exporting VST values
vst_values_mcfkd <- vst(mcfkd.de,blind=FALSE)
vst_values_scckd <- vst(scckd.de,blind=FALSE)
vst_values_koscc <- vst(koscc.de,blind=FALSE)
write.csv(assay(vst_values_mcfkd), file="deseq2_vst_MCF10A_p63_knockdowns.csv")
write.csv(assay(vst_values_scckd), file="deseq2_vst_SCC25_p63_knockdowns.csv")
write.csv(assay(vst_values_koscc), file="deseq2_vst_MCF10Ap53KO_vs_SCC25.csv")
#get DEG results table comparing different cell types
res.MCF10KD <- results(mcfkd.de, contrast=c("genotype", "MCF10A_p63sh","MCF10A_SCR"))
res.SCC25KD <- results(scckd.de, contrast=c("genotype", "SCC25_p63_shRNA","SCC25_ctr_shRNA"))
res.MCF10Ap53KOvsSCC25 <- results(koscc.de, contrast=c("genotype", "SCC25_ctr_shRNA","MCF10Ap53KO_DMSO"))
# create a list of ensembl genes in the v104 database
# add gene names and symbols in the results lists to make downstream analysis easier
ensembl_symbols <- ensembl_gene_list %>%
                  dplyr::select(ensembl_gene_id,hgnc_symbol,chromosome_name,start_position,
                                end_position,strand,transcript_is_canonical) %>%
                  distinct() %>% filter(transcript_is_canonical == 1)
#convert DEseq2 results tables to tibbles
tb.MCF10Ap53KOvsSCC25 = res.MCF10Ap53KOvsSCC25 %>% data.frame() %>%
                        rownames_to_column(var="ensembl_gene_id") %>% as_tibble()
tb.SCC25KD = res.SCC25KD %>% data.frame() %>% rownames_to_column(var="ensembl_gene_id") %>% as_tibble()
tb.MCF10KD = res.MCF10KD %>% data.frame() %>% rownames_to_column(var="ensembl_gene_id") %>% as_tibble()
# using left_join to merge DEseq2 results tables with ENSEMBL geneID and hgnc symbol information
tb.MCF10Ap53KOvsSCC25 <- left_join(tb.MCF10Ap53KOvsSCC25,ensembl_symbols,by="ensembl_gene_id")
tb.SCC25KD <- left_join(tb.SCC25KD,ensembl_symbols,by="ensembl_gene_id")
tb.MCF10KD = left_join(tb.MCF10KD,ensembl_symbols,by="ensembl_gene_id")
#write to csv
write.csv(tb.MCF10Ap53KOvsSCC25, 'deseq2_MCF10Ap53KO_DMSO_vs_SCC25_ctrl_shRNA.csv')
write.csv(tb.SCC25KD, 'deseq2_SCC25_shRNA_knockdown.csv')
write.csv(tb.MCF10KD, 'deseq2_MCF10A_shRNA_knockdown.csv')