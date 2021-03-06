### Genomic Features of Muscle-Invasive Bladder Cancer Arising After Prostate Radiotherapy
## December 2021

library(tidyverse)
library(dplyr)
library(SigMA)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(NMF)
library(survival)
library(survminer)
library(ggpubr)
library(wesanderson)


## Survival analysis Fig 1
# OS Mouw and TCGA 
mouw_surv <- read.csv("survival_data.csv", header = TRUE)

mouw_surv$dead_binary <- as.numeric(mouw_surv$dead_binary)
mouw_surv$time_os..days. <- as.numeric(mouw_surv$time_os..days.)
mouw_surv$time_prog..days. <- as.numeric(mouw_surv$time_prog..days.)
mouw_surv$recurrence.progression <- as.numeric(mouw_surv$recurrence.progression)

mouw_surv <- dplyr::rename(mouw_surv, type= tumor_sample_barcode, OS = dead_binary, OS.time= time_os..days., 
                           DFI= recurrence.progression, DFI.time= time_prog..days.) %>%
        dplyr::select(type, OS, OS.time, DFI, DFI.time)


tcga_clinic_info <- read.csv("TCGA_CDR.csv", header = TRUE)
tcga_bca_clinic <- filter(tcga_clinic_info, type == "BLCA") %>%
        dplyr::select(bcr_patient_barcode, type, OS, OS.time, DSS, DSS.time, DFI, DFI.time)

surv_tcga <- dplyr::select(tcga_bca_clinic, type, OS, OS.time, DFI, DFI.time) 
surv_tcga$OS <- as.numeric(surv_tcga$OS)
surv_tcga$OS.time <- as.numeric(surv_tcga$OS.time)
surv_tcga$DFI <- as.numeric(surv_tcga$DFI)
surv_tcga$DFI.time <- as.numeric(surv_tcga$DFI.time)

surv_tcga_mouw <- bind_rows(surv_tcga, mouw_surv)  

surv_object <- Surv(time = surv_tcga_mouw$OS.time, event = surv_tcga_mouw$OS)
surv_object
os_surv <- survfit(surv_object ~ type, data = surv_tcga_mouw)
summary(os_surv)
plotOS <- ggsurvplot(os_surv, data = surv_tcga_mouw, pval = TRUE, risk.table = TRUE, tables.height = 0.2,
                     title= "", font.main =20,
                     font.x = 20,
                     font.y = 20,
                     ylab= "OS",
                     xlab= "Days from surgery",
                     font.tickslab = 20,
                     break.time.by= 1000,
                     legend.labs=c("TCGA", "RA BC"), 
                     legend.title='', legend = c(0.8,0.8), palette = wes_palette("FantasticFox1", n=2), 
                     fontsize = 10,
                     tables.theme = clean_theme(),
                     font.title.risk.table= 10,
                     risk.table.y.text = T, 
                     risk.table.font.y  = 200
)
plotOS$plot <- plotOS$plot +
        theme(legend.text = element_text(size = 20))
plotOS

## =============================================================================
## Analysis of variants, indels and copy-number alteration
# Path to Mouw Radiation-associated BCa MAF files

input_path <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/"
input_suffix <- "*.common_variant_filter.pass.maf"

file_list_passmaf <- list.files(path = input_path, pattern = input_suffix)
file_list_passmaf <- paste(input_path, file_list_passmaf, sep = "")

file_list_passmaf

df=data.frame()
for(i in 1:length(file_list_passmaf)) {
        df_tmp <- read.table(file_list_passmaf[i], sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))
        df <- rbind(df,df_tmp)
        rm(df_tmp)
}

df <- mutate(df, total_reads = df$t_alt_count + df$t_ref_count)
df_reduced <- subset(df, total_reads >= 14)


## Define common mutated genes in bladder cancer and DNA repair genes -- Comutation plots generated with Comut tool
bca_genes <- c("TP53", "KMT2D", "KDM6A", 'ARID1A',"PIK3CA", "KMT2C", "RB1", "EP300", "FGFR3", "STAG2", "ATM", "FAT1", "ELF3", "CREBBP", "ERBB2", "SPTAN1", "KMT2A", "ERBB3", "ERCC2", "CDKN1A", "ASXL2", "TSC1", "FBXW7") 

dnarepair_genes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "BRCA1", "BRCA2", "MRE11A", "NBN", "RAD50", "RAD51", "RAD51B", "RAD51D", "RAD52", "RAD54L", "BRIP1", "FANCA", "FANCC", "PALB2", "RAD51C", "BLM", "ATM", "ATR", "CHEK1", "CHEK2", "MDC1", "POLE", "MUTYH", "PARP1", "RECQL4")

## Copy-number analysis
# Create a dataframes of seg files from all samples
tumor_seg=data.frame()
for(i in 1:length(file_list_tumorseg)) {
        df_tmp <- read.table(file_list_tumorseg[i], sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))
        tumor_seg <- rbind(tumor_seg,df_tmp)
        rm(df_tmp)
}

# Rename the column names
colnames(tumor_seg) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Probes", "Segment_Mean")

# Export the dataframe as a tab txt file for GISTIC
write.table(tumor_seg, file= "tumor_seg_22samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)
 
# Plot CDKNA in comut (Fig. 2A)
cna_mouw <- read.csv("all_lesions_99.txt", sep= "\t", head= TRUE)
cna_cdkn2a_mouw <- cna_mouw %>%
        dplyr::slice(39) %>%
        t()
write.table(cna_cdkn2a_mouw, file= "bladder_cdkn2a.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# Genome plot (Suppl Fig. 5) generated in maftools
mouw.gistic <- readGistic(gisticAllLesionsFile = "all_lesions_99.txt", gisticAmpGenesFile = "amp_genes_99.txt", gisticDelGenesFile = "del_genes_99.txt", gisticScoresFile = "tumor_gistic_scores.gistic", isTCGA = FALSE)
gisticChromPlot(gistic = mouw.gistic, markBands = "all", ref.build = "hg19", cytobandOffset = 0.1, 
                txtSize = 0.8, cytobandTxtSize = 1)

## =============================================================================
## Mutational Signatures  -- Fig 3A
# Sigma
maf.file <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/df_reduced.tsv"
maf <- fread(maf.file)
maf$Chromosome <- as.character(maf$Chromosome)
maf[which(maf$Chromosome == "23"), ]$Chromosome <- "X"
maf[which(maf$Chromosome == "24"), ]$Chromosome <- "Y"
write.table(maf, "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/df_reduced.tsv", row.names = F, sep = '\t', quote = F)
maf.file <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/df_reduced.tsv"
genomes_matrix <- SigMA:::make_matrix(maf.file, file_type = 'maf', ref_genome_name = "hg19")
genomes <- SigMA:::conv_snv_matrix_to_df(genomes_matrix)
write.table(genomes,
            "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/df_reduced.tsv",
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)
SigMA:::run("/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/df_reduced.tsv",
            "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/SigMA_w_MSI_output.txt",
            data = "seqcap", 
            do_assign = T,
            tumor_type = "bladder",
            do_mva = T, 
            lite_format = T, 
            check_msi = T)

Mouw_SigMA <- read.csv(file= "~/VAL/Mouw_RA_MIBC/eur_urol/rt_maf/SigMA/SigMA_w_MSI_output.txt", header = TRUE, sep = ",")
str(Mouw_SigMA)
write.table(Mouw_SigMA, file= "Mouw_SigMA_october.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Plot signatures by SNV in each sample
Mouw_SigMA$exps_all
r <- strsplit(Mouw_SigMA$exps_all, split = "[_]")
df_test = data.frame(tumor = rep(Mouw_SigMA$tumor, sapply(r, length)), exps_all = unlist(r))

s <- strsplit(Mouw_SigMA$sigs_all, split = "[.]")
df_plot = data.frame(tumor = rep(Mouw_SigMA$tumor, sapply(s, length)), sigs_all = unlist(s)) %>%
        dplyr::rename(Tumor = tumor)

# Combine two tables above and round SNVs
sigs_snv <- cbind(df_test, df_plot) %>%
        dplyr::select(tumor, sigs_all, exps_all) 
sigs_snv$exps_all <- as.numeric(sigs_snv$exps_all) 
sigs_snv <- sigs_snv %>% mutate_at(3, funs(round(., 0)))
write.table(sigs_snv, file= "sigs_snv_october.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## Change Signature_1 by Sig1 etc
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_1")] <- "Sig1"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_2")] <- "Sig2"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_3")] <- "Sig3"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_5")] <- "Sig5"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_13")] <- "Sig13"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_18")] <- "Sig18"

# Reorder signatures for boxplot
sigs_snv$sigs_all <- factor(sigs_snv$sigs_all, levels = c("Sig1", "Sig2", "Sig13", "Sig3", "Sig5", "Sig18"))

# Remove samples with sig18 and plot
no_sig18 <- sigs_snv %>% filter(sigs_all != "Sig18")
no_sig18$sigs_all <- factor(no_sig18$sigs_all, levels = c("Sig1", "Sig2", "Sig13", "Sig3", "Sig5"))
no_sig18
plotno_sig18 <- no_sig18 %>% ggplot(aes(x= sigs_all, y = exps_all)) +
        geom_boxplot(aes(fill=sigs_all)) +
        scale_fill_manual(values = c("yellow", "orange","gray", "brown", "#56B4E9", "#999999")) +
        xlab("Mutational Signatures") +
        ylab("Number of SNVs") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20, angle = 90), axis.title.x = element_text(colour="black", size = 20, margin = margin(t=20, r=0, b=0, l=0)), axis.title.y = element_text(colour="black", size = 20, margin = margin(t=0, r=20, b=0, l=0))) 

plotno_sig18 

# Pairwise comparison between Signatures
compare_means(exps_all ~ sigs_all, data = no_sig18, paired = FALSE, p.adjust.method = "bonferroni")

## =============================================================================
## To make appropriate comparisions between radiation-associated and non-radiation-associated bladder cancer cohorts, we created an intersection interval and bait list, and proceed both cohorts through the tumor only pipeline

# Radiation-associated MIBC
mouw_xrt_i <- read.table("xrt-tumor-only_common-variant-filter.txt", sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA")) # need to remove the rows that equal 1 for the column common_variant
mouw_xrt <- filter(mouw_xrt_i, common_variant == '0')
mouw_xrt <- mutate(mouw_xrt, total_reads = mouw_xrt$t_alt_count + mouw_xrt$t_ref_count)
mouw_xrt <- subset(mouw_xrt, total_reads >= 14)

# Non-radiation associated from tumor only pipeline
ercc2_to_i <- read.table("ercc2-tumor-only_common-variant-filter.txt", sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))  # need to remove the rows that equal 1 for the column common_variant
ercc2_to <- filter(ercc2_to_i, common_variant == '0')
ercc2_to <- mutate(ercc2_to, total_reads = ercc2_to$t_alt_count + ercc2_to$t_ref_count)
ercc2_to <- subset(ercc2_to, total_reads >= 14)


# Remove noncoding variants and keep Missense, Nonsense, Splice Site, FrameShift indels, indels and nonsop mutations
mouw_xrt <- mouw_xrt %>% filter(Variant_Classification != "Intron", Variant_Classification != "Silent", Variant_Classification != "intron", Variant_Classification != "3'UTR", Variant_Classification != "5'UTR", Variant_Classification != "IGR", Variant_Classification != "RNA", Variant_Classification != "5'Flank", Variant_Classification != "5'utr", Variant_Classification != "Start_Codon_Del", Variant_Classification != "igr", Variant_Classification != "lincRNA", Variant_Classification != "3'utr", Variant_Classification != "Start_Codon_SNP", Variant_Classification != "De_novo_Start_InFrame", Variant_Classification != "De_novo_Start_OutOfFrame", Variant_Classification != "5'flank") 

# Same for non-radiated tumors
ercc2_to <- ercc2_to %>% filter(Variant_Classification != "Intron", Variant_Classification != "Silent", Variant_Classification != "intron", Variant_Classification != "3'UTR", Variant_Classification != "5'UTR", Variant_Classification != "IGR", Variant_Classification != "RNA", Variant_Classification != "5'Flank", Variant_Classification != "5'utr", Variant_Classification != "Start_Codon_Del", Variant_Classification != "igr", Variant_Classification != "lincRNA", Variant_Classification != "3'utr", Variant_Classification != "Start_Codon_SNP", Variant_Classification != "De_novo_Start_InFrame", Variant_Classification != "De_novo_Start_OutOfFrame", Variant_Classification != "5'flank") 


variants_tmb <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation","Translation_Start_Site",
                  "Nonstop_Mutation")
variants_snv <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation","Translation_Start_Site",
                  "Nonstop_Mutation", "Silent")
variants_allindel <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Stop_Codon_Del", "Start_Codon_Del", "Start_Codon_Ins", "Stop_Codon_Ins")
variants_allinser <- c("Frame_Shift_Ins", "In_Frame_Ins", "Stop_Codon_Ins", "Start_Codon_Ins")
variants_alldel <- c("Frame_Shift_Del", "In_Frame_Del","Stop_Codon_Del", "Start_Codon_Del")


## Fig 2A left pannel -- TMB
tmb_xrt <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)

tmb_ercc2_to <-  ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)

## for boxplots
tmb_ercc2_to$case.id[c(1:50)] <- "Non-radiated"
tmb_xrt$case.id[c(1:22)] <- "Radiation-associated"

## TMB for tumor only cohorts
tmb_to_cohorts <- bind_rows(tmb_xrt, tmb_ercc2_to)
tmb_to_cohorts$case.id <- factor(tmb_to_cohorts$case.id, levels = c("Non-radiated", "Radiation-associated"))

# plot deletion per megabase 
plot_tmb_to <- tmb_to_cohorts %>% ggplot(aes(x= case.id, y = tmb)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("TMB / Mb") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  

plot_tmb_to

wilcox.test(tmb~case.id, data = tmb_to_cohorts) 

#### Fig 3C left panel -- ratio indels/ SNVs

allindel_xrt <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_allindel)%>%
        filter(`Variant_Classification %in% variants_allindel`== "TRUE")%>%
        select(Tumor_Sample_Barcode, allindel = n)

snv_xrt_t <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_snv)%>%
        filter(`Variant_Classification %in% variants_snv`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, allsnv = n)

allindel_ercc2_to <- ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_allindel)%>%
        filter(`Variant_Classification %in% variants_allindel`== "TRUE")%>%
        select(Tumor_Sample_Barcode, allindel = n)

snv_ercc2_to_t <- ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_snv)%>%
        filter(`Variant_Classification %in% variants_snv`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, allsnv = n)

# for boxplot indels/snv 
xrt_indel_variants_t <- data.frame(left_join(allindel_xrt, snv_xrt_t, by= "Tumor_Sample_Barcode") %>% 
                                           mutate(ratio = allindel/allsnv))
xrt_indel_variants_t <- dplyr::rename(xrt_indel_variants_t, case.id = Tumor_Sample_Barcode)

ercc2_indel_variants_t <- data.frame(left_join(allindel_ercc2_to, snv_ercc2_to_t, by= "Tumor_Sample_Barcode") %>% 
                                             mutate(ratio = allindel/allsnv))
ercc2_indel_variants_t <- dplyr::rename(ercc2_indel_variants_t,case.id = Tumor_Sample_Barcode)

ercc2_indel_variants_t$case.id[c(1:50)] <- "Non-radiated"
xrt_indel_variants_t$case.id[c(1:22)] <- "Radiation-associated"

indel_snvs_t <- bind_rows(xrt_indel_variants_t, ercc2_indel_variants_t)
indel_snvs_t$case.id <- factor(indel_snvs_t$case.id, levels = c("Non-radiated", "Radiation-associated"))

# Plots indels /snv
plotindelsnv_t <- indel_snvs_t %>% ggplot(aes(x= case.id, y = ratio)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Ratio Indels / SNVs") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  
plotindelsnv_t

wilcox.test(ratio~case.id, data = indel_snvs_t) 


#### Fig 3C middle panel --- Del/SNV
snv_xrt_t <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_snv)%>%
        filter(`Variant_Classification %in% variants_snv`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, allsnv = n)

alldel_xrt <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n)

snv_ercc2_to_t <- ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_snv)%>%
        filter(`Variant_Classification %in% variants_snv`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, allsnv = n)

alldel_ercc2_to <- ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n)

# Compute ratios del/snv 
xrt_del_variants_t <- data.frame(left_join(alldel_xrt, snv_xrt_t, by= "Tumor_Sample_Barcode") %>% 
                                         mutate(ratio = alldel/allsnv))
xrt_del_variants_t <- dplyr::rename(xrt_del_variants_t, case.id = Tumor_Sample_Barcode)

ercc2_del_variants_t <- data.frame(left_join(alldel_ercc2_to, snv_ercc2_to_t, by= "Tumor_Sample_Barcode") %>% 
                                           mutate(ratio = alldel/allsnv))
ercc2_del_variants_t <- dplyr::rename(ercc2_del_variants_t, case.id = Tumor_Sample_Barcode)

# for boxplots
ercc2_del_variants_t$case.id[c(1:50)] <- "Non-radiated"
xrt_del_variants_t$case.id[c(1:22)] <- "Radiation-associated"

ratios_dels_snv3_t <- bind_rows(xrt_del_variants_t, ercc2_del_variants_t)
ratios_dels_snv3_t$case.id <- factor(ratios_dels_snv3_t$case.id, levels = c("Non-radiated", "Radiation-associated"))

# Plots del / snv 
plotdelsnv3_t <- ratios_dels_snv3_t %>% ggplot(aes(x= case.id, y = ratio)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Ratio Deletions / SNVs") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  
plotdelsnv3_t

wilcox.test(ratio~case.id, data = ratios_dels_snv3_t) 

#### Fig 3C right panel -- deletions per Mb 
xrt_del_mb <- mouw_xrt%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        mutate(del_mb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n, del_mb)
xrt_del_mb <- dplyr::rename(xrt_del_mb, case.id = Tumor_Sample_Barcode) 

ercc2_to_del_mb <- ercc2_to%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        mutate(del_mb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n, del_mb)
ercc2_to_del_mb <- dplyr::rename(ercc2_to_del_mb, case.id = Tumor_Sample_Barcode) 

ercc2_to_del_mb$case.id[c(1:50)] <- "Non-radiated"
xrt_del_mb$case.id[c(1:22)] <- "Radiation-associated"

ratios_del_mb_to <- bind_rows(xrt_del_mb, ercc2_to_del_mb)
ratios_del_mb_to$case.id <- factor(ratios_del_mb_to$case.id, levels = c("Non-radiated", "Radiation-associated"))

# Plot deletions per megabase
plotdelmb_to <- ratios_del_mb_to %>% ggplot(aes(x= case.id, y = del_mb)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Deletions per Mb") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  
plotdelmb_to 

wilcox.test(del_mb~case.id, data = ratios_del_mb_to) 

## Plots for fig3 C 
ggarrange(plotindelsnv_t, NULL, plotdelsnv3_t, NULL, plotdelmb_to,
          nrow = 1, widths = c(1, 0.75, 1, 0.75, 1), 
          legend = "right", common.legend = TRUE )


## =============================================================================
#### Suppl fig 1 -- Survival analysis Radiation-associated cohort        
# for OS          
surv_object <- Surv(time = mouw_surv$OS.time, event = mouw_surv$OS)
surv_object
mouw_os_surv <- survfit(surv_object ~ type, data = mouw_surv)
summary(mouw_os_surv)
plotmow_os <- ggsurvplot(mouw_os_surv, data = mouw_surv, pval = TRUE, risk.table = TRUE, tables.height = 0.2,
                         conf.int = FALSE,
                         legend = "none",
                         title= "", font.main =14,
                         font.x = 20,
                         font.y = 20,
                         ylab= "OS",
                         xlab= "Days from surgery",
                         font.tickslab = 20,
                         break.time.by= 1000,
                         fontsize = 10,
                         font.title.risk.table= 8, 
                         font.risk.table= 8, 
                         risk.table.y.text = FALSE, risk.table.x.text = FALSE,
)
plotmow_os
survfit(Surv(mouw_surv$OS.time, mouw_surv$OS) ~ 1, data = mouw_surv)

# for RFS
rfs_object <- Surv(time = mouw_surv$DFI.time, event = mouw_surv$DFI)
rfs_object
rfs_surv <- survfit(rfs_object ~ type, data = mouw_surv)
summary(rfs_surv)
plotmouw_rfs <- ggsurvplot(rfs_surv, data = mouw_surv, pval = TRUE, risk.table = TRUE, tables.height = 0.2,
                           conf.int = FALSE,
                           legend = "none",
                           title= "", font.main =14,
                           font.x = 20,
                           font.y = 20,
                           ylab= "RFS",
                           xlab= "Days from surgery",
                           font.tickslab = 20,
                           break.time.by= 1000,
                           fontsize = 10,
                           font.title.risk.table= 10,
                           font.risk.table= 10, 
                           risk.table.y.text = FALSE, risk.table.x.text = FALSE,
)
plotmouw_rfs

## =============================================================================
#### Suppl fig 3 and 7
## Plots for with fraction of patients with mutations in commonly mutated genes in MIBC. 
mut_mibcgenes_mouw <- mouw_xrt%>%
        dplyr::filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

mut_mibcgenes_ercc2 <- ercc2_to %>%
        dplyr::filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

## Example for TP53
dplyr::count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("TP53"), ])
dplyr::count(mut_mibcgenes_ercc2[mut_mibcgenes_ercc2$Hugo_Symbol== c("TP53"), ]) 
fraction_TP53 <- data.frame(
        Cohort = c("Non-Radiated", "Non-Radiated", "RA BC", "RA BC"),
        Fraction = c((50-28)/50, 28/50, ((19-12)/19), 12/19),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_TP53 <- matrix(c((50-28), 28, ((19-12)), 12), nrow = 2)
fisher.test(to_fisher_TP53) # p-value = 0.7855
plottp53 <- ggplot(data = fraction_TP53, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("TP53") +
        ylab("Fraction of patients") + xlab("Cohort") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 14), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 14), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 14))
