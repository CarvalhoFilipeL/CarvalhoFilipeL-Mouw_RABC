#### Genomic Features of Muscle-Invasive Bladder Tumors Arising After Prostate Radiotherapy
#### June 2021

## Packages
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(wesanderson)
library(maftools)
library(survival)
library(survminer)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(NMF)
library(SigMA)
library(data.table)

## Survival analysis figure 1
setwd("/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/km_tables/")
mouw_surv <- read.csv("survival_data.csv", header = TRUE)

mouw_surv$dead_binary <- as.numeric(mouw_surv$dead_binary)
mouw_surv$time_os..days. <- as.numeric(mouw_surv$time_os..days.)
mouw_surv$time_prog..days. <- as.numeric(mouw_surv$time_prog..days.)
mouw_surv$recurrence.progression <- as.numeric(mouw_surv$recurrence.progression)

mouw_surv <- dplyr::rename(mouw_surv, type= tumor_sample_barcode, OS = dead_binary, OS.time= time_os..days., 
                           DFI= recurrence.progression, DFI.time= time_prog..days.) %>%
        dplyr::select(type, OS, OS.time, DFI, DFI.time)

# mouw cohort alone          
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

# OS Mouw and TCGA (Suppl Fig 1) 
setwd("/Users/filipecarvalho/VAL/TCGA/")
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
                     title= "", font.main =14,
                     font.x = 14,
                     font.y = 14,
                     ylab= "OS",
                     xlab= "Days from surgery",
                     font.tickslab = 12,
                     break.time.by= 1000,
                     legend.labs=c("TCGA", "RA BC"), 
                     legend.title="", legend = c(0.8,0.83), palette = c("darkblue", "brown"),
                     fontsize = 4,
                     font.title.risk.table= 4 )
plotOS

## ==================================================================================
## Analyses of variantes and indels
setwd("/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/maf/")

input_path <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/maf/"
input_suffix <- "*.commonfilter.pass.maf"

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
write.table(df_reduced, file= "df_reduced.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## ==================================================================================
## Copy-number analysis
## Create a dataframes of seg files from all samples
tumor_seg=data.frame()
for(i in 1:length(file_list_tumorseg)) {
        df_tmp <- read.table(file_list_tumorseg[i], sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))
        tumor_seg <- rbind(tumor_seg,df_tmp)
        rm(df_tmp)
}

## Rename the column names
colnames(tumor_seg) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Probes", "Segment_Mean")

## Export the dataframe as a tab txt file for GISTIC
write.table(tumor_seg, file= "tumor_seg_22samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## Processing copy-number data from GISTIC in maftools
## reading and summarizing gistic output files
mouw.gistic <- readGistic(gisticAllLesionsFile = "all_lesions_99.txt", gisticAmpGenesFile = "amp_genes_99.txt", gisticDelGenesFile = "del_genes_99.txt", gisticScoresFile = "tumor_gistic_scores.gistic", isTCGA = FALSE)

## GISTIC object
mouw.gistic ## Samples 22, nGenes 3653, cytoBands  45, Amp  50323, Del    847, total 51170   

## Plot CDKNA in comut (Fig. 2)
cna_mouw <- read.csv("all_lesions_99.txt", sep= "\t", head= TRUE)
cna_cdkn2a_mouw <- cna_mouw %>%
        dplyr::slice(39) %>%
        t()
write.table(cna_cdkn2a_mouw, file= "bladder_cdkn2a.csv", sep = "\t", row.names = FALSE, quote = FALSE)

## Genome plot (Suppl Fig. 6)
gisticChromPlot(gistic = mouw.gistic, markBands = "all", ref.build = "hg19", cytobandOffset = 0.1, 
                txtSize = 0.8, cytobandTxtSize = 1)

## ==================================================================================
## Mutational signatures (Fig. 3A)
maf.file <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/maf/df_reduced.tsv"
maf <- fread(maf.file)
maf$Chromosome <- as.character(maf$Chromosome)
maf[which(maf$Chromosome == "23"), ]$Chromosome <- "X"
maf[which(maf$Chromosome == "24"), ]$Chromosome <- "Y"
write.table(maf, "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/df_reduced.tsv", row.names = F, sep = '\t', quote = F)
maf.file <- "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/df_reduced.tsv"
genomes_matrix <- SigMA:::make_matrix(maf.file, file_type = 'maf', ref_genome_name = "hg19")
genomes <- SigMA:::conv_snv_matrix_to_df(genomes_matrix)
write.table(genomes,
            "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/df_reduced.tsv",
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)
SigMA:::run("/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/df_reduced.tsv",
            "/Users/filipecarvalho/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/SigMA_w_MSI_output.txt",
            data = "seqcap", 
            do_assign = T,
            tumor_type = "bladder",
            do_mva = T, 
            lite_format = T, 
            check_msi = T)

Mouw_SigMA <- read.csv(file= "~/VAL/Mouw_RA_MIBC/Mouw_WES_RA_MIBC.Rproj/data/SigMA/SigMA_w_MSI_output.txt", header = TRUE, sep = ",")
str(Mouw_SigMA)

### Plot signatures by SNV in each sample
Mouw_SigMA$exps_all
r <- strsplit(Mouw_SigMA$exps_all, split = "[_]")
df_test = data.frame(tumor = rep(Mouw_SigMA$tumor, sapply(r, length)), exps_all = unlist(r))

s <- strsplit(Mouw_SigMA$sigs_all, split = "[.]")
df_plot = data.frame(tumor = rep(Mouw_SigMA$tumor, sapply(s, length)), sigs_all = unlist(s)) %>%
        dplyr::rename(Tumor = tumor)

# combine two tables above and round SNVs
sigs_snv <- cbind(df_test, df_plot) %>%
        dplyr::select(tumor, sigs_all, exps_all) 
sigs_snv$exps_all <- as.numeric(sigs_snv$exps_all) 
sigs_snv <- sigs_snv %>% mutate_at(3, funs(round(., 0)))

## change Signature_1 by Sig1 etc
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_1")] <- "Sig1"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_2")] <- "Sig2"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_3")] <- "Sig3"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_5")] <- "Sig5"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_13")] <- "Sig13"
sigs_snv$sigs_all[which(sigs_snv$sigs_all == "Signature_18")] <- "Sig18"

# reorder signatures for boxplot
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
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20, angle = 90), axis.title.x = element_text(colour="black", size = 20, margin = margin(t=20, r=0, b=0, l=0)), axis.title.y = element_text(colour="black", size = 20, margin = margin(t=0, r=20, b=0, l=0))) +
        stat_compare_means(method = "kruskal.test")
plotno_sig18 

# pairwise comparison between Signatures
compare_means(exps_all ~ sigs_all, data = no_sig18, paired = FALSE, p.adjust.method = "bonferroni")

## ==================================================================================
## ## Analyses SNVs and indels in RA-MIBC and TCGA (Fig. 3C, 3D, and Supplementary Fig. 4)
## read subset of bladder cancer mafs from TCGA mc3 
setwd("/Users/filipecarvalho/VAL/TCGA/")

tcga_bca_maf <- read.table("mc3.blca.maf", header = TRUE, sep = "")

# Compute and compare genes mutated in primary MIBC and DNA mismatch repair genes in TCGA vs RA-MIBC (Fig 3C and Supplementary Fig. 4)
bca_genes <- c("TP53", "KMT2D", "KDM6A", 'ARID1A',"PIK3CA", "KMT2C", "RB1", "EP300", "FGFR3", "STAG2", "ATM", "FAT1", "ELF3", "CREBBP", "ERBB2", "SPTAN1", "KMT2A", "ERBB3", "ERCC2", "CDKN1A", "ASXL2", "TSC1", "FBXW7") 

dnarepair_genes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "BRCA1", "BRCA2", "MRE11A", "NBN", "RAD50", "RAD51", "RAD51B", "RAD51D", "RAD52", "RAD54L", "BRIP1", "FANCA", "FANCC", "PALB2", "RAD51C", "BLM", "ATM", "ATR", "CHEK1", "CHEK2", "MDC1", "POLE", "MUTYH", "PARP1", "RECQL4")

mut_mibcgenes_mouw <- df_reduced%>%
        dplyr::filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

mut_mibcgenes_tcga <- tcga_bca_maf %>%
        dplyr::filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

mut_drrgenes_mouw <- df_reduced%>%
        dplyr::filter(Hugo_Symbol%in%dnarepair_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

mut_drrgenes_tcga <- tcga_bca_maf %>%
        dplyr::filter(Hugo_Symbol%in%dnarepair_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        dplyr::distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)

# Examples of barplots in Fig 3C and Suppl Fig 4. Remaining plots were generated with similar code.
# FANCA
dplyr::count(mut_drrgenes_mouw[mut_drrgenes_mouw$Hugo_Symbol== c("FANCA"), ]) #5
dplyr::count(mut_drrgenes_tcga[mut_drrgenes_tcga$Hugo_Symbol== c("FANCA"), ]) #30
fraction_FANCA <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-30)/411, 30/411, ((19-5)/19), 5/19),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_FANCA <- matrix(c((411-30), 30, ((19-5)), 5), nrow = 2)
fisher.test(to_fisher_FANCA) # p-value = 0.01349
plotfanca <- ggplot(data = fraction_FANCA, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("FANCA") +
        ylab("Fraction of patients") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 16), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 16), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 16))
plotfanca
# TP53
dplyr::count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("TP53"), ])#12
dplyr::count(mut_mibcgenes_tcga[mut_mibcgenes_tcga$Hugo_Symbol== c("TP53"), ]) #203
fraction_TP53 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-203)/411, 203/411, ((19-12)/19), 12/19),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_TP53 <- matrix(c((411-203), 203, ((19-12)), 12), nrow = 2)
fisher.test(to_fisher_TP53) # p-value = 0.3483
plottp53 <- ggplot(data = fraction_TP53, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("TP53") +
        ylab("Fraction of patients") + xlab("Cohort") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 14), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 14), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 14))
plottp53

## Analyses of SNVs and indels in RA-MIBC as imprint of expouse to radiation therapy (Fig 3D and Suppl Fig 3)
variants_all <- unique(df_reduced$Variant_Classification)
variants_tmb <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation","Translation_Start_Site",
                  "Nonstop_Mutation")
variants_alldel <- c("Frame_Shift_Del", "In_Frame_Del","Stop_Codon_Del", "Start_Codon_Del")

tmb_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)
snv_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_all)%>%
        filter(`Variant_Classification %in% variants_all`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, allsnv = n)
alldel_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n)

tmb_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(patient.id, non_silent_muts = n, tmb)
snv_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_all)%>%
        filter(`Variant_Classification %in% variants_all`== "TRUE")%>%
        dplyr::select(patient.id, allsnv = n)
alldel_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        dplyr::select(patient.id, alldel = n)

#### Compute and compare deletions per Mb in RA-MIBC and TCGA
mouw_del_mb <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        mutate(del_mb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, alldel = n, del_mb)
mouw_del_mb <- dplyr::rename(mouw_del_mb, case.id = Tumor_Sample_Barcode) 

tcga_del_mb <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        mutate(del_mb= n/30)%>%
        dplyr::select(patient.id, alldel = n, del_mb)
tcga_del_mb <- dplyr::rename(tcga_del_mb, case.id = patient.id)

## boxplot deletions per Mb
tcga_del_mb$case.id[c(1:367)] <- "TCGA"
mouw_del_mb$case.id[c(1:22)] <- "RA BC"
ratios_del_mb <- bind_rows(mouw_del_mb, tcga_del_mb)
ratios_del_mb$case.id <- factor(ratios_del_mb$case.id, levels = c("TCGA", "RA BC"))
plotdelmb <- ratios_del_mb %>% ggplot(aes(x= case.id, y = del_mb)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Deletions per Mb") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# Stats with Wilcoxon rank sum/Mann Whitney
wilcox.test(del_mb~case.id, data = ratios_del_mb) #  p-value = 9.394e-14

## Compute and compare ratios Del/SNVs in RA-MIBC and TCGA 
mouw_del_variants <- data.frame(left_join(alldel_mouw, snv_mouw, by= "Tumor_Sample_Barcode") %>% 
                                        mutate(ratio = alldel/allsnv))
mouw_del_variants <- dplyr::rename(mouw_del_variants, case.id = Tumor_Sample_Barcode)

tcga_del_variants <- data.frame(left_join(alldel_tcga, snv_tcga, by= "patient.id") %>% 
                                        mutate(ratio = alldel/allsnv))
tcga_del_variants <- dplyr::rename(tcga_del_variants,case.id = patient.id)

## boxplots Del/SNVs
tcga_del_variants$case.id[c(1:367)] <- "TCGA"
mouw_del_variants$case.id[c(1:22)] <- "RA BC"
ratios_dels_snv <- bind_rows(mouw_del_variants, tcga_del_variants)
ratios_dels_snv$case.id <- factor(ratios_dels_snv$case.id, levels = c("TCGA", "RA BC"))

# Plots del / snv
plotdelsnv <- ratios_dels_snv %>% ggplot(aes(x= case.id, y = ratio)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Ratio Deletions / SNVs") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

## Stats with Wilcoxon rank sum/Mann Whitney
wilcox.test(ratio~case.id, data = ratios_dels_snv) # p-value = 2.475e-07

ggarrange(plotdelmb, NULL, plotdelsnv,
          nrow = 1, widths = c(1, 0.25, 1), 
          legend = "right", common.legend = TRUE 
          )

## Calculate and plot TMB in RA-MIBC and TCGA (Suppl. Fig 3)
tmb_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)

tmb_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(patient.id, non_silent_muts = n, tmb)

##for TMB boxplots
tmb_tcga <- dplyr::rename(tmb_tcga,case.id = patient.id)
tmb_tcga$case.id[c(1:411)] <- "TCGA"
tmb_mouw <- dplyr::rename(tmb_mouw,case.id = Tumor_Sample_Barcode)
tmb_mouw$case.id[c(1:22)] <- "RA BC"
#plot_tmb_tcga <- dplyr::select(tmb_tcga, case.id, tmb)
#plot_tmb_mouw <- dplyr::select(tmb_mouw, case.id, tmb)
plot_tmb <- bind_rows(plot_tmb_tcga, plot_tmb_mouw)
plot_tmb$case.id <- factor(plot_tmb$case.id, levels = c("TCGA", "RA BC"))

# Plots TMB 
plottmb <- plot_tmb %>% ggplot(aes(x= case.id, y = tmb)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("TMB / Mb") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))
plottmb

## Stats with Wilcoxon rank sum/Mann Whitney
wilcox.test(tmb~case.id, data = plot_tmb) # W = 2513, p-value = 0.0004472

## ==================================================================================
## CNA Genome plot Suppl Fig 6
# Reading and summarizing gistic output files
mouw.gistic <- readGistic(gisticAllLesionsFile = "all_lesions_99.txt", gisticAmpGenesFile = "amp_genes_99.txt", gisticDelGenesFile = "del_genes_99.txt", gisticScoresFile = "tumor_gistic_scores.gistic", isTCGA = FALSE)

# GISTIC object
mouw.gistic ## Samples 22, nGenes 3653, cytoBands  45, Amp  50323, Del    847, total 51170   

# Genome plot
gisticChromPlot(gistic = mouw.gistic, markBands = "all", ref.build = "hg19", cytobandOffset = 0.1, 
                txtSize = 0.8, cytobandTxtSize = 1)