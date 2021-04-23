### Mutational Features of Muscle-Invasive Bladder Tumors Arising After Prostate Radiotherapy
### WES sequencing, Illumina, hg19, Tumor-only, n=27, FFPE samples. 

# Packages
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


## Survival analysis figure 1
mouw_surv <- read.csv("survival_data.csv", header = TRUE)

mouw_surv$dead_binary <- as.numeric(mouw_surv$dead_binary)
mouw_surv$time_os..days. <- as.numeric(mouw_surv$time_os..days.)
mouw_surv$time_prog..days. <- as.numeric(mouw_surv$time_prog..days.)
mouw_surv$recurrence.progression <- as.numeric(mouw_surv$recurrence.progression)
mouw_surv <- mouw_surv %>% 
        select(-(mrn))
mouw_surv <- dplyr::rename(mouw_surv, type= tumor_sample_barcode, OS = dead_binary, OS.time= time_os..days., 
                           DFI= recurrence.progression, DFI.time= time_prog..days.) %>%
        select(type, OS, OS.time, DFI, DFI.time)

# combine OS and DFI mouw and TCGA
str(surv_tcga)

surv_tcga <- select(tcga_bca_clinic, type, OS, OS.time, DFI, DFI.time) 
surv_tcga$OS <- as.numeric(surv_tcga$OS)
surv_tcga$OS.time <- as.numeric(surv_tcga$OS.time)
surv_tcga$DFI <- as.numeric(surv_tcga$DFI)
surv_tcga$DFI.time <- as.numeric(surv_tcga$DFI.time)

surv_tcga_mouw <- bind_rows(surv_tcga, mouw_surv)  
str(surv_tcga_mouw)

# surv TCGA and Mouw
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
# for RFS
rfs_object <- Surv(time = surv_tcga_mouw$DFI.time, event = surv_tcga_mouw$DFI)
rfs_object
rfs_surv <- survfit(rfs_object ~ type, data = surv_tcga_mouw)
summary(rfs_surv)
plotrfs <- ggsurvplot(rfs_surv, data = surv_tcga_mouw, pval = TRUE, risk.table = TRUE, tables.height = 0.2,
                      title= "", font.main =14,
                      font.x = 14,
                      font.y = 14,
                      ylab= "RFS",
                      xlab= "Days from surgery",
                      font.tickslab = 12,
                      break.time.by= 1000,
                      legend.labs=c("TCGA", "RA BC"), 
                      legend.title="", legend = c(0.8,0.83), palette = c("darkblue", "brown"),
                      fontsize = 4,
                      font.title.risk.table= 4 )

plotrfs

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

## ==============================================================================
## Analyses of variants and indels

setwd("/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/rt_pass_mafs/")

input_path <- "/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/rt_pass_mafs/"
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
dim(df_reduced) 

## =============================================================================
## Mutational signatures 
xrt_bca_maf <- read.maf(maf = df_reduced, verbose = FALSE)
xrt_bca_tnm <- trinucleotideMatrix(maf = xrt_bca_maf, prefix = "chr", add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

xrt_bca_sig <- estimateSignatures(mat = xrt_bca_tnm, nTry = 5)
plotCophenetic(res = xrt_bca_sig) 

## cophenetic plot showed biggest drop n=3
xrt_bca_sig <- extractSignatures(mat = xrt_bca_tnm, n= 3)

## Compare against COSMIC version3 60 signatures
xrt_bca_v3_cosm <- compareSignatures(nmfRes = xrt_bca_sig, sig_db = "SBS")

## Plots COSMIC signatures
pheatmap::pheatmap(mat = xrt_bca_v3_cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

maftools::plotSignatures(nmfRes = xrt_bca_sig, title_size = 1.5, sig_db = "SBS")

## ===================================================================================
## Maftools df_reduced for ALL DNA damage response genes
dnarepair_genes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "BRCA1", "BRCA2", "MRE11A", "NBN", "RAD50", "RAD51", "RAD51B", "RAD51D", "RAD52", "RAD54L", "BRIP1", "FANCA", "FANCC", "PALB2", "RAD51C", "BLM", "ATM", "ATR", "CHEK1", "CHEK2", "MDC1", "POLE", "MUTYH", "PARP1", "RECQL4")

## RColorBrewer color pallet for variant classifications
# vc_col <- RColorBrewer::brewer.pal(n=6, name = "Dark2")
# names(vc_col) <- c("Frame_Shift_Ins", "Frame_Shift_Del", "Missense_Mutation", "Splice_Site", "Nonsense_Mutation", "Multi_Hit")
# print(vc_col)
# oncoplot(maf = xrt_bca_maf,  colors= vc_col, genes = dnarepair_genes, draw_titv = TRUE)

## mutated DNA repair genes in MOUW ----> FANCA, MSH6, BRACA1, CHEK2, RAD54L, PALB2, POLE, ERCC5, ATM, MSH2, ERCC2, BRCA2, RAD51B, MDC1

### ============================================================================
#read subset of bladder cancer mafs from TCGA mc3 and clinical info
tcga_bca_maf <- read.table("mc3.blca.maf", header = TRUE, sep = "")
tcga_clinic_info <- read.csv("TCGA_CDR.csv", header = TRUE)

# str(tcga_bca_maf)
# tcga_clinic_info$type

tcga_bca_all <- filter(tcga_clinic_info, type == "BLCA")
write.table(tcga_bca_all, file= "tcga_bca_all.csv", sep = "\t", row.names = FALSE, quote = FALSE)

tcga_bca_clinic <- filter(tcga_clinic_info, type == "BLCA") %>%
                        select(bcr_patient_barcode, type, OS, OS.time, DSS, DSS.time, DFI, DFI.time)

## Analyses SNV/indels in RA-MIBC and TCGA
variants_all <- unique(df_reduced$Variant_Classification)
variants_tmb <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation","Translation_Start_Site",
                  "Nonstop_Mutation")
variants_truncating <- c("Frame_Shift_Del", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Ins", "Translation_Start_Site","Nonstop_Mutation")
variants_indel <- c("In_Frame_Del", "In_Frame_Ins")
variants_silent <- c("3'UTR",  "Silent", "Intron", "5'UTR", "RNA", "5'Flank", "3'Flank")

variants_allindel <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Stop_Codon_Del", "Start_Codon_Del", "Start_Codon_Ins", "Stop_Codon_Ins")
variants_fsindel <- c("Frame_Shift_Del", "Frame_Shift_Ins")
variants_allinser <- c("Frame_Shift_Ins", "In_Frame_Ins", "Stop_Codon_Ins", "Start_Codon_Ins")
variants_alldel <- c("Frame_Shift_Del", "In_Frame_Del","Stop_Codon_Del", "Start_Codon_Del")

bca_genes <- c("TP53", "KMT2D", "KDM6A", "PIK3CA", "KMT2C", "RB1", "EP300", "FGFR3", "STAG2", "ATM", "FAT1", "ELF3", "CREBBP", "ERBB2", "SPTAN1", "KMT2A", "ERBB3", "ERCC2", "CDKN1A", "ASXL2", "TSC1", "FBXW7")
dnarepair_genes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "BRCA1", "BRCA2", "MRE11A", "NBN", "RAD50", "RAD51", "RAD51B", "RAD51D", "RAD52", "RAD54L", "BRIP1", "FANCA", "FANCC", "PALB2", "RAD51C", "BLM", "ATM", "ATR", "CHEK1", "CHEK2", "MDC1", "POLE", "MUTYH", "PARP1", "RECQL4")

# number of variants
mouw_silent <- df_reduced%>%
        dplyr::count(Variant_Classification %in% variants_silent)
mouw_missense <- df_reduced%>%
        dplyr::count(Variant_Classification == 'Missense_Mutation')
mouw_nonsense <- df_reduced%>%
        dplyr::count(Variant_Classification == 'Nonsense_Mutation')
mouw_nonsense <- df_reduced%>%
        dplyr::count(Variant_Classification %in% variants_fsindel)

#Compute and compare TMB and all indels, SNV/dels, del/ins ratios
tmb_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)
write.table(tmb_mouw, file= "tmb_mouw.csv", sep = "\t", row.names = FALSE, quote = FALSE)

allindel_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_allindel)%>%
        filter(`Variant_Classification %in% variants_allindel`== "TRUE")%>%
        select(Tumor_Sample_Barcode, allindel = n)
snv_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_all)%>%
        filter(`Variant_Classification %in% variants_all`== "TRUE")%>%
        select(Tumor_Sample_Barcode, allsnv = n)
alldel_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        select(Tumor_Sample_Barcode, alldel = n)
allins_mouw <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_allinser)%>%
        filter(`Variant_Classification %in% variants_allinser`== "TRUE")%>%
        select(Tumor_Sample_Barcode, allinser = n)

## TCGA
tmb_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        select(patient.id, non_silent_muts = n, tmb)
snv_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_all)%>%
        filter(`Variant_Classification %in% variants_all`== "TRUE")%>%
        select(patient.id, allsnv = n)
allindel_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_allindel)%>%
        filter(`Variant_Classification %in% variants_allindel`== "TRUE")%>%
        select(patient.id, allindel = n)

alldel_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_alldel)%>%
        filter(`Variant_Classification %in% variants_alldel`== "TRUE")%>%
        select(patient.id, alldel = n)
allins_tcga <- tcga_bca_maf%>%
        group_by(patient.id)%>%
        dplyr::count(Variant_Classification %in% variants_allinser)%>%
        filter(`Variant_Classification %in% variants_allinser`== "TRUE")%>%
        select(patient.id, allinser = n)

## Compute ratios indels/snv 
mouw_indel_variants <- data.frame(left_join(allindel_mouw, snv_mouw, by= "Tumor_Sample_Barcode") %>% 
                       mutate(ratio = allindel/allsnv))
mouw_indel_variants <- dplyr::rename(mouw_indel_variants, case.id = Tumor_Sample_Barcode)

tcga_indel_variants <- data.frame(left_join(allindel_tcga, snv_tcga, by= "patient.id") %>% 
                        mutate(ratio = allindel/allsnv))
tcga_indel_variants <- dplyr::rename(tcga_indel_variants,case.id = patient.id)

## for boxplots
tcga_indel_variants$case.id[c(1:386)] <- "TCGA"
mouw_indel_variants$case.id[c(1:27)] <- "RA BC"

ratio_tcga <- select(tcga_indel_variants, case.id, ratio)
ratio_mouw <- select(mouw_indel_variants, case.id, ratio)
ratios <- bind_rows(ratio_mouw, ratio_tcga)
ratios$case.id <- factor(ratios$case.id, levels = c("TCGA", "RA BC"))
levels(ratios$case.id)
ratios

# Plots indels /snv
plotsnvratio <- ratios %>% ggplot(aes(x= case.id, y = ratio)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Ratio Indels / SNVs") +
        theme_bw() + theme(legend.title = element_blank(),  legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))
plotsnvratio

## Stats with Wilcoxon rank sum/Mann Whitney
wilcox.test(ratio~case.id, data = ratios) # W = 9048, p-value = 1.573e-10

## Compute ratios del/ins
mouw_del_inser <- data.frame(left_join(alldel_mouw, allins_mouw, by= "Tumor_Sample_Barcode") %>% 
                                     mutate(ratio.indel = alldel/allinser))
mouw_del_inser <- dplyr::rename(mouw_del_inser, case.id = Tumor_Sample_Barcode) 

tcga_del_inser <- data.frame(left_join(alldel_tcga, allins_tcga, by= "patient.id") %>% 
                                    drop_na(alldel, allinser) %>%
                                      mutate(ratio.indel = alldel/allinser))
tcga_del_inser <- dplyr::rename(tcga_del_inser, case.id = patient.id)

## for boxplots
tcga_del_inser$case.id[c(1:278)] <- "TCGA"
mouw_del_inser$case.id[c(1:27)] <- "RA BC"

ratio_indels_tcga <- select(tcga_del_inser, case.id, ratio.indel)
ratio_indels_mouw <- select(mouw_del_inser, case.id, ratio.indel)
ratios_indels <- bind_rows(ratio_indels_mouw, ratio_indels_tcga)

ratios_indels$case.id <- factor(ratios_indels$case.id, levels = c("TCGA", "RA BC"))
levels(ratios_indels$case.id)

# Plots del / insertion
plotindels <- ratios_indels %>% ggplot(aes(x= case.id, y = ratio.indel)) +
        geom_boxplot(aes(fill=case.id)) +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        xlab("Cohorts") +
        ylab("Ratio Deletions / Insertions") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))
plotindels

## Stats with Wilcoxon rank sum/Mann Whitney
wilcox.test(ratio.indel~case.id, data = ratios_indels) #W = 4295.5, p-value = 0.2136

##for TMB boxplots
tmb_tcga <- dplyr::rename(tmb_tcga,case.id = patient.id)
tmb_tcga$case.id[c(1:411)] <- "TCGA"

tmb_mouw <- dplyr::rename(tmb_mouw,case.id = Tumor_Sample_Barcode)
tmb_mouw$case.id[c(1:27)] <- "RA BC"

plot_tmb_tcga <- select(tmb_tcga, case.id, tmb)
plot_tmb_mouw <- select(tmb_mouw, case.id, tmb)
plot_tmb <- bind_rows(plot_tmb_tcga, plot_tmb_mouw)

#sanity check
plot_tmb$case.id <- factor(plot_tmb$case.id, levels = c("TCGA", "RA BC"))
levels(plot_tmb$case.id)
plot_tmb

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
wilcox.test(tmb~case.id, data = plot_tmb) # W = 2888, p-value = 2.982e-05

## ==============================================================================
## plot MIBC mutated genes in TCGA vs Mouw
mut_mibcgenes_mouw <- df_reduced%>%
        filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)
dim(mut_mibcgenes_mouw) #99 2
mut_mibcgenes_mouw
mut_mibcgenes_tcga <- tcga_bca_maf %>%
        filter(Hugo_Symbol%in%bca_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)
dim(mut_mibcgenes_tcga) # 1543  2

unique(tcga_bca_maf$Tumor_Sample_Barcode) #411
unique(df_reduced$Tumor_Sample_Barcode) # 27

## plot mutated DNA repair genes in TCGA vs Mouw. Given presence Sig6 in RA MIBC cohort, are there higher frequency of mutation in Mouw vs TCGA?
mut_drrgenes_mouw <- df_reduced%>%
        filter(Hugo_Symbol%in%dnarepair_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)
dim(mut_drrgenes_mouw) #65  2

mut_drrgenes_tcga <- tcga_bca_maf %>%
        filter(Hugo_Symbol%in%dnarepair_genes)%>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        distinct() %>%
        dplyr::rename(case.id = Tumor_Sample_Barcode)
dim(mut_drrgenes_tcga) # 581   2

# FANCA
count(mut_drrgenes_mouw[mut_drrgenes_mouw$Hugo_Symbol== c("FANCA"), ]) #8
count(mut_drrgenes_tcga[mut_drrgenes_tcga$Hugo_Symbol== c("FANCA"), ]) #30
fraction_FANCA <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-30)/411, 30/411, ((27-8)/27), 8/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_FANCA <- matrix(c((411-30), 30, ((27-8)), 8), nrow = 2)
fisher.test(to_fisher_FANCA) # p-value = 0.001007

plotfanca <- ggplot(data = fraction_FANCA, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("FANCA") +
        ylab("Fraction of patients") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# BRCA1
count(mut_drrgenes_mouw[mut_drrgenes_mouw$Hugo_Symbol== c("BRCA1"), ]) # 5
count(mut_drrgenes_tcga[mut_drrgenes_tcga$Hugo_Symbol== c("BRCA1"), ]) # 22
fraction_BRCA1 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-22)/411, 22/411, ((27-5)/27), 5/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_BRCA1 <- matrix(c((411-22), 11, ((27-5)), 5), nrow = 2)
fisher.test(to_fisher_BRCA1) # p-value = 0.001869

plotbrca1 <- ggplot(data = fraction_BRCA1, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("BRCA1") +
        ylab("Fraction of patients") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 16), legend.title = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 16), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 16))

# MSH6
count(mut_drrgenes_mouw[mut_drrgenes_mouw$Hugo_Symbol== c("MSH6"), ]) #6
count(mut_drrgenes_tcga[mut_drrgenes_tcga$Hugo_Symbol== c("MSH6"), ]) #14
fraction_MSH6 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-14)/411, 14/411, ((27-6)/27), 6/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_MSH6 <- matrix(c((411-14), 14, ((27-6)), 6), nrow = 2)
fisher.test(to_fisher_MSH6) # p-value = 0.0006645

plotmsh6 <- ggplot(data = fraction_MSH6, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("MSH6") +
        ylab("Fraction of patients") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# CHEK2
count(mut_drrgenes_mouw[mut_drrgenes_mouw$Hugo_Symbol== c("CHEK2"), ]) #5
count(mut_drrgenes_tcga[mut_drrgenes_tcga$Hugo_Symbol== c("CHEK2"), ]) #11
fraction_CHEK2 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-11)/411, 11/411, ((27-5)/27), 5/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_CHEK2 <- matrix(c((411-11), 11, ((27-5)), 5), nrow = 2)
fisher.test(to_fisher_CHEK2) # p-value = 0.001666

plotchek2 <- ggplot(data = fraction_CHEK2, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("CHEK2") +
        ylab("Fraction of patients") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# Figure with all plots together for fig 3C
ggarrange(plotfanca, plotbrca1, 
          plotmsh6, plotchek2,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = 'right')

## ========================================================================================
## Examples of barplots in the suppl fig 5. Remaining plots were generated with similar code
# TP53
count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("TP53"), ]) #14
count(mut_mibcgenes_tcga[mut_mibcgenes_tcga$Hugo_Symbol== c("TP53"), ]) #203
fraction_TP53 <- data.frame(
                Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
                Fraction = c((411-203)/411, 203/411, ((27-14)/27), 14/27),
                Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_TP53 <- matrix(c((411-203), 203, ((27-14)), 14), nrow = 2)
fisher.test(to_fisher_TP53) # p-value = 0.8446

plottp53 <- ggplot(data = fraction_TP53, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
                geom_bar(stat= "identity") +
                ggtitle("TP53") +
                ylab("Fraction of patients") + xlab("Cohort") +
                scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
                scale_y_continuous(expand = c(0,0)) +
                theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# RB1
count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("RB1"), ]) #6
count(mut_mibcgenes_tcga[mut_mibcgenes_tcga$Hugo_Symbol== c("RB1"), ]) #78
fraction_RB1 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-78)/411, 78/411, ((27-6)/27), 6/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_RB1 <- matrix(c((411-78), 78, ((27-6)), 6), nrow = 2)
fisher.test(to_fisher_RB1) # p-value = 0.6205

plotrb1 <- ggplot(data = fraction_RB1, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("RB1") +
        ylab("Fraction of patients") + xlab("Cohort") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# ERBB2
count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("ERBB2"), ]) #5
count(mut_mibcgenes_tcga[mut_mibcgenes_tcga$Hugo_Symbol== c("ERBB2"), ]) #61
fraction_ERBB2 <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-61)/411, 61/411, ((27-5)/27), 5/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_ERBB2 <- matrix(c((411-61), 61, ((27-5)), 5), nrow = 2)
fisher.test(to_fisher_ERBB2) # p-value = 0.5805

ploterbb2 <- ggplot(data = fraction_ERBB2, mapping = aes(x = Cohort, y = Fraction, fill=Mutation))+
        geom_bar(stat= "identity") +
        ggtitle("ERBB2") +
        ylab("Fraction of patients") + xlab("Cohort") +
        scale_fill_manual(values = wes_palette(n =2, name = "FantasticFox1")) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20), legend.title = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20))

# for KDM6A
count(mut_mibcgenes_mouw[mut_mibcgenes_mouw$Hugo_Symbol== c("KDM6A"), ]) #11
count(mut_mibcgenes_tcga[mut_mibcgenes_tcga$Hugo_Symbol== c("KDM6A"), ]) #113
fraction_KDM6A <- data.frame(
        Cohort = c("TCGA", "TCGA", "RA BC", "RA BC"),
        Fraction = c((411-113)/411, 113/411, ((27-11)/27), 11/27),
        Mutation = c("Wt", "Mut", "Wt", "Mut"))
to_fisher_KDM6A <- matrix(c((411-113), 113, ((27-11)), 11), nrow = 2)
fisher.test(to_fisher_KDM6A) 


## ==============================================================================
## Calculate TMB mutated samples in Mouw
tmb_mouw_silent_synonymous <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_silent)%>%
        filter(`Variant_Classification %in% variants_silent`== "TRUE")  %>%
        select(Tumor_Sample_Barcode, silent_muts = n) %>%
        dplyr::rename(Synonymous = silent_muts)

tmb_mouw_nonsynonymous <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        select(Tumor_Sample_Barcode, non_silent_muts = n) %>%
        dplyr::rename(Nonsynonymous = non_silent_muts)

mouw_silent_nonsil <- merge(tmb_mouw_nonsynonymous, tmb_mouw_silent_synonymous, by= "Tumor_Sample_Barcode")
write.table(mouw_silent_nonsil, file= "mouw_sile_nonsil.csv", sep = "\t", row.names = FALSE, quote = FALSE)

## ==============================================================================
## CNA
## Seg files from "cheungatm/GATK4_Somatic_CNV_workflow_w_plots" pipeline
setwd("/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/seg_27samples/")

input_path <- "/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/seg_27samples/"
input_suffix <- "*.seg"

file_list_tumorseg <- list.files(path = input_path, pattern = input_suffix)
file_list_tumorseg <- paste(input_path, file_list_tumorseg, sep = "")

### Create a dataframes for tumor_seg_29samples
tumor_seg=data.frame()
for(i in 1:length(file_list_tumorseg)) {
        df_tmp <- read.table(file_list_tumorseg[i], sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))
        tumor_seg <- rbind(tumor_seg,df_tmp)
        rm(df_tmp)
}

## Rename the column names
colnames(tumor_seg) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Probes", "Segment_Mean")

## Export the dataframe as a tab txt file for GISTIC
write.table(tumor_seg, file= "tumor_seg_27samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Set working directory for GISTIC output files
setwd("/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/seg_27samples/")

## Processing copy-number data from GISTIC in maftools
## reading and summarizing gistic output files
mouw.gistic <- readGistic(gisticAllLesionsFile = "all_lesions_99.txt", gisticAmpGenesFile = "amp_genes_99.txt", gisticDelGenesFile = "del_genes_99.txt", gisticScoresFile = "tumor_gistic_scores.gistic", isTCGA = FALSE)

## GISTIC object
mouw.gistic ## Samples 27, nGenes  3567, cytoBands  49, Amp   65209, Del    705, total  65914 

## Genome plot
gisticChromPlot(gistic = mouw.gistic, markBands = "all", ref.build = "hg19", cytobandOffset = 0.1, 
                txtSize = 0.8, cytobandTxtSize = 1)

## Plot CDKNA in comut
cna_mouw <- read.csv("gistic_all_lesions.txt", sep= "\t", head= TRUE)
cna_cdkn2a_mouw <- cna_mouw %>%
        dplyr::slice(42)
write.table(cna_cdkn2a_mouw, file= "bladder_cdkn2a.csv", sep = "\t", row.names = FALSE, quote = FALSE)

### ============================================================================
## Paired samples. Create df of paired samples for comut
# unique(df_reduced[c('Tumor_Sample_Barcode')])
pair_samples <- c("RP-2203_BS10A12568-1_v3_Exome_OnPrem", "RP-2203_BS09N20349-1_v3_Exome_OnPrem",
                  "RP-2203_BS12A07361-1_v3_Exome_OnPrem", "RP-2203_BS11M23658-1_v3_Exome_OnPrem",
                  "RP-2203_BS17J06357-1_v3_Exome_OnPrem", "RP-2203_BS16X63420-1_v3_Exome_OnPrem", "RP-2203_BS17R22951-1_v3_Exome_OnPrem")
df_pairs <- df_reduced%>%
        filter(Tumor_Sample_Barcode %in% pair_samples) 

write.table(df_pairs, file= "df_pairs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### ============================================================================
## Supplementary tables gene list of variants and indels
str(df_reduced)

mut_bca_genes <- c('ASXL2','CDKN1A','ERCC2','CREBBP','PIK3CA', 'KMT2A', 'ATM','KMT2C','KMT2D','ERBB2','STAG2','SPTAN1','ELF3','EP300','RB1','FGFR3','FBXW7','TSC1','FAT1', 'ERBB3','KDM6A','TP53')

install.packages("writexl")
library(writexl)
write_xlsx(suppl_bcagene_list,"/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/suppl_bcagenes_list.xlsx")

suppl_bcagene_list <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        select(Tumor_Sample_Barcode, Hugo_Symbol, Entrez_Gene_Id, Chromosome, Start_position, End_position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Transcript_Strand, cDNA_Change, Protein_Change, t_alt_count, t_ref_count) %>%
        dplyr::filter(Hugo_Symbol %in% mut_bca_genes)

write_xlsx(suppl_mmrgene_list,"/Users/filipecarvalho/RT_MIBC_Mouw_2020/data/suppl_mmrgene_list.xlsx")

suppl_mmrgene_list <- df_reduced%>%
        group_by(Tumor_Sample_Barcode)%>%
        select(Tumor_Sample_Barcode, Hugo_Symbol, Entrez_Gene_Id, Chromosome, Start_position, End_position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Transcript_Strand, cDNA_Change, Protein_Change, t_alt_count, t_ref_count) %>%
        dplyr::filter(Hugo_Symbol %in% dnarepair_genes)