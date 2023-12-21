library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(gridExtra)
library(vegan)
library(corrplot)
library(RColorBrewer)
library(Hmisc)
library(ggExtra)
library(FactoMineR)
library(factoextra)
library(cluster)
library(gplots)
library(rstatix)
library(pheatmap)
library(ggrepel)
library(sva)
library(santaR)
library(tempted)
library(impute)
library(ggVennDiagram)
library(mixOmics)
library(ggpubr)

setwd('C:/Users/sthomas/OneDrive - University of California, San Diego Health/Mouse Antibiotics/All')

# Cleanup feature quant table from MZmine ###########################################################################
# norm <- read.csv("Mouse_Antibiotics_iimn_GNPS_quant.csv")
# norm <- norm %>% dplyr::select(-X)
# norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ""))
# norm <- norm %>% dplyr::select(-(row.m.z:correlation.group.ID), -(best.ion:neutral.M.mass))
# norm <- norm %>% mutate(across(row.ID:annotation.network.number, as.character))
# 
# ## Merge ions from the same molecule (iimn)
# iimn <- norm %>%
#   group_by(annotation.network.number) %>%
#   summarise(across(where(is.numeric), ~sum(., na.rm = TRUE))) %>%
#   filter(!is.na(annotation.network.number)) %>%
#   rename(row.ID = annotation.network.number)
# labels <- norm %>% dplyr::select(row.ID, annotation.network.number)
# iimn$row.ID <- paste0(iimn$row.ID, "_i")
# norm_iimn <- norm %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
# norm_iimn <- rbind(norm_iimn, iimn)
# 
# ## Read in FBMN library IDs from GNPS
# library_ID <- read.delim("FBMN_All.tsv")
# library_ID <- library_ID %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex, RTConsensus)
# colnames(library_ID) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number", "RTConsensus")
# library_ID$Precursor_Mass <- round(library_ID$Precursor_Mass, 2)
# library_ID$row.ID <- as.character(library_ID$row.ID)
# library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")
# library_ID[library_ID == "N/A"] <- NA
# library_ID$Network_Number[library_ID$Network_Number == -1] <- NA
# 
# # Standard Test
# test <- full_join(library_ID, norm, by = "row.ID") 
# test <- test %>% dplyr::select(row.ID:RTConsensus, contains("CA"))
# 
# ## Collapse ions (keeps all IDs, precursor masses, and network numbers of features in the iimn network separated by OR)
# result_list <- list()
# for (colname in c("LibraryID", "Precursor_Mass", "Network_Number")) {
#   result <- library_ID %>%
#     filter(!is.na(annotation.network.number)) %>%
#     group_by(annotation.network.number, !!sym(colname)) %>%
#     summarise(n = n()) %>%
#     group_by(annotation.network.number) %>%
#     summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>%
#     rename(!!colname := "concat")
#   result_list[[colname]] <- result
# }
# ## Calculate mean RT for all iimn ions
# for (colname in c("RTConsensus")) {
#   result <- library_ID %>%
#     filter(!is.na(annotation.network.number)) %>%
#     group_by(annotation.network.number) %>%
#     summarise(n = mean(!!sym(colname))) %>%
#     rename(!!colname := "n")
#   result_list[[colname]] <- result
# }
# 
# library_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
# library_iimn$row.ID <- paste0(library_iimn$row.ID, "_i")
# library_f <- library_ID %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
# library_iimn <- rbind(library_f, library_iimn)
# library_iimn <- library_iimn %>% rename(Metabolite = row.ID)
# library_iimn[library_iimn == ""] <- NA
# write_csv(library_iimn, "library_matches_all.csv")
# 
# ## Do the same for the bile acid library
# library_ID <- read.delim("FBMN_Bile.tsv")
# library_ID <- library_ID %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex, RTConsensus)
# colnames(library_ID) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number", "RTConsensus")
# library_ID$Precursor_Mass <- round(library_ID$Precursor_Mass, 2)
# library_ID$row.ID <- as.character(library_ID$row.ID)
# library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")
# library_ID[library_ID == "N/A"] <- NA
# library_ID$Network_Number[library_ID$Network_Number == -1] <- NA
# 
# # Look for standards 
# library_test <- read.delim("FBMN_Bile_Acyl.tsv")
# library_test <- library_test %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex, RTConsensus)
# colnames(library_test) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number", "RTConsensus")
# library_test$Precursor_Mass <- round(library_test$Precursor_Mass, 2)
# library_test$row.ID <- as.character(library_test$row.ID)
# library_test[library_test == "N/A"] <- NA
# test <- full_join(library_test, norm, by = "row.ID") 
# test <- test %>% dplyr::select(row.ID:RTConsensus, contains("CA"))
# 
# ## Collapse ions (keeps all IDs, precursor masses, and network numbers of features in the iimn network separated by OR)
# result_list <- list()
# for (colname in c("LibraryID", "Precursor_Mass", "Network_Number")) {
#   result <- library_ID %>%
#     filter(!is.na(annotation.network.number)) %>%
#     group_by(annotation.network.number, !!sym(colname)) %>%
#     summarise(n = n()) %>%
#     group_by(annotation.network.number) %>%
#     summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>%
#     rename(!!colname := "concat")
#   result_list[[colname]] <- result
# }
# ## Calculate mean RT for all iimn ions
# for (colname in c("RTConsensus")) {
#   result <- library_ID %>%
#     filter(!is.na(annotation.network.number)) %>%
#     group_by(annotation.network.number) %>%
#     summarise(n = mean(!!sym(colname))) %>%
#     rename(!!colname := "n")
#   result_list[[colname]] <- result
# }
# 
# library_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
# library_iimn$row.ID <- paste0(library_iimn$row.ID, "_i")
# library_f <- library_ID %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
# library_iimn <- rbind(library_f, library_iimn)
# library_iimn <- library_iimn %>% rename(Metabolite = row.ID)
# library_iimn[library_iimn == ""] <- NA
# write_csv(library_iimn, "library_matches_bile.csv")
# 
# norm_iimn <- norm_iimn %>%
#   gather(Sample, value, 2:ncol(norm_iimn)) %>%
#   spread(row.ID, value)
# norm_iimn$Sample <- str_replace_all(norm_iimn$Sample, "\\.", "-")
# norm_iimn <- column_to_rownames(norm_iimn, var = "Sample")
# ## Filter metabolites that are found in at least 10% of samples
# norm_iimn <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.1*nrow(norm_iimn))))
# norm_cleaned <- norm_iimn %>% rownames_to_column("Sample")
# ## Plot TICs and geometric mean of all samples (this can help identify any poor quality runs)
# TICs <- norm_cleaned
# TICs[TICs == 0] <- NA
# geom_mean <- function(x, na.rm = FALSE) {
#   x <- x[x > 0]  # Exclude non-positive values
#   if (length(x) > 0) exp(mean(log(x), na.rm = na.rm)) else NA
# }
# TICs <- TICs %>% rowwise() %>% mutate(TIC = sum(c_across(where(is.numeric)), na.rm = TRUE), Mean = geom_mean(c_across(where(is.numeric)), na.rm = TRUE))
# TIC_plot <- TICs %>% dplyr::select(Sample, TIC, Mean) %>% pivot_longer(TIC:Mean, names_to = "Type", values_to = "value")
# ggplot(TIC_plot, aes(x=reorder(Sample, value), y=log2(value), color=Type)) +
#   geom_point() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# TICs[is.na(TICs)] <- 0
# ## Create dataframe with features normalized by total ion current (TIC)
# TICs <- TICs %>%  mutate(across(where(is.numeric), ~ ./TIC)) %>% dplyr::select(-TIC, -Mean)
# ## Write csv of un-preprocessed data and TIC-normalized data
# write_csv(norm_cleaned, "Metabolites_cleaned.csv")
# write_csv(TICs, "Metabolites_TIC.csv")
# 
# ## Perform rclr preprocessing
# norm_iimn <- norm_iimn %>% decostand(method = "rclr")
# norm_iimn <- norm_iimn %>% rename_with(~str_replace_all(., "^", "X"))
# 
# ## Missing value imputation
# ## For this dataframe, we'll only take metabolites found in 20% of samples
# norm_imp <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.2*nrow(norm_iimn))))
# norm_imp <- norm_imp %>% t() %>% as.data.frame()
# ## We also need to filter out any samples with more than 80% missing data
# norm_imp <- norm_imp %>% dplyr::select(where(~sum(. != 0) >= (0.25*nrow(norm_imp))))
# norm_imp[norm_imp == 0] <- NA
# norm_imp <- impute.knn(as.matrix(norm_imp), rowmax=1, k=50, rng.seed=1234)
# norm_imp <- norm_imp$data %>% t() %>% as.data.frame()
# 
# norm_iimn <- rownames_to_column(norm_iimn, var = "Sample")
# norm_imp <- rownames_to_column(norm_imp, var = "Sample")
# 
# write_csv(norm_iimn, "Metabolites_normalized.csv")
# write_csv(norm_imp, "Metabolites_normalized_imputed.csv")
# write_csv(labels, "iimn_groups.csv")
 
# ## Combat to look at batch effects ##############################################################
# metadata <- read.delim("Metadata_All.tsv")
# metadata <- metadata %>% mutate(across(2:ncol(.), as.factor))
# norm <- read.csv("Metabolites_normalized_imputed.csv")
# batch <- metadata %>% filter(Type == "Fecal")
# batch$Batch <- paste0(batch$Hypothesis, "_", batch$Batch)
# data <- norm %>%
#   gather(Metabolite, value, 2:ncol(norm)) %>%
#   spread(Sample, value)
# data <- column_to_rownames(data, var = "Metabolite")
# batch <- batch %>% filter(Sample %in% colnames(data))
# data <- data %>% dplyr::select(any_of(batch$Sample))
# mod <- model.matrix(~AMP+Mother_Infant, data=batch)
# 
# data_corr <- ComBat(dat = data, batch = batch$Batch, mod=mod, par.prior = TRUE)
# data_corr <- as.data.frame(data_corr)
# data_corr <- rownames_to_column(data_corr, var = "Metabolite")
# data <- data_corr %>%
#   gather(Sample, value, 2:ncol(data_corr)) %>%
#   spread(Metabolite, value)
# write_csv(data, "Metabolites_normalized_imputed_bc.csv")
# 
# ## Distance matrix to check how well it worked ######################################################
# norm <- read.csv("Metabolites_TIC.csv")
# data <- inner_join(metadata, norm, by = "Sample")
# data1 <- data %>%
#   filter(Type == "Fecal") %>% filter(Mother_Infant == "Mother") %>%
#   arrange(Hypothesis)
# data1 <- column_to_rownames(data1, var = "Sample")
# data2 <- data1 %>%
#   dplyr::select(-(Hypothesis:Experiment))
# 
# # Plot distance matrix. This will show how similar each of your sample sets are to each other. 0 = most similar
# res.dist <- get_dist(data2, method = "spearman")
# fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"), lab_size = 5, order = FALSE)
# 
# # Plot PCA by mass spec run to look for batch effects
# pca <- PCA(data1[,11:ncol(data1)], scale.unit = FALSE, ncp = 5, graph = FALSE)
# fviz_screeplot(pca)
# 
# # Now you can plot the first two principle components
# fviz_pca_ind(pca, axes = c(1,2), habillage = as.factor(data1$Run), label = "ind", addEllipses = TRUE, ellipse.level = 0.95, invisible="quali")


## You only need to do normalization once ######################################################
## If you've already done it, you can skip to this step 
metadata <- read.delim("Metadata_All.tsv")
metadata$host_subject_id <- paste0(metadata$Cage, "_", metadata$Hypothesis, "_", metadata$Batch)
metadata$Experiment <- paste0(metadata$Hypothesis, "_", metadata$AMP)
bile_lib <- read.csv("library_matches_bile.csv")
bile_lib <- bile_lib %>% filter(!is.na(LibraryID)) %>% dplyr::select(-Network_Number)
bile_lib$Metabolite <- paste0("X", bile_lib$Metabolite)
bile_lib$Bile_Type <- ifelse(str_detect(bile_lib$LibraryID, "Mono"), "Mono", (ifelse(str_detect(bile_lib$LibraryID, "Di"), "Di", (ifelse(str_detect(bile_lib$LibraryID, "Tri"), "Tri", 
                            (ifelse(str_detect(bile_lib$LibraryID, "Non"), "Non", (ifelse(str_detect(bile_lib$LibraryID, "Tetra"), "Tetra", "Penta")))))))))
bile_lib$Bile_Type <- factor(bile_lib$Bile_Type, levels = c("Non", "Mono", "Di", "Tri", "Tetra", "Penta"))

norm <- read.csv("Metabolites_normalized.csv")
norm_imp <- read.csv("Metabolites_normalized_imputed_bc.csv")
cleaned <- read.csv("Metabolites_cleaned.csv")
TIC <- read.csv("Metabolites_TIC.csv")
data <- inner_join(metadata, norm, by = "Sample") 
data <- data %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite))

## Create venn diagram with amount of bile acids and blood and feces
venn <- data %>% pivot_longer(cols =c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
venn$Expression[venn$Expression == 0] <- NA
venn <- venn %>%  filter(!is.na(Expression))
venn <- inner_join(venn, bile_lib, by = c("Bile_Acids" = "Metabolite"), multiple = "any")

## Big fig
test <- read.csv("Metabolites_cleaned.csv")
test <- inner_join(metadata, test, by = "Sample")
test <- test %>% pivot_longer(cols =c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
test$Expression[test$Expression == 0] <- NA
test <- test %>%  filter(!is.na(Expression))
test <- inner_join(test, bile_lib, by = c("Bile_Acids" = "Metabolite"), multiple = "any")
test$LibraryID <- str_remove_all(test$LibraryID, " OR .*")
test$LibraryID <- str_remove_all(test$LibraryID, "Candidate ")
test$LibraryID <- str_remove_all(test$LibraryID, "sirius.*")
test$LibraryID <- str_remove_all(test$LibraryID, "buddy.*")
test <- test %>% filter(Bile_Acids == "X2051_i")
ggplot(test, aes(x=Days_Post_Birth, y=Expression)) +
  geom_point(aes(color = as.factor(AMP))) +
  #theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~LibraryID, scales = "free_y") 

## Serum vs. Fecal
x <- list(A = venn %>% filter(Type == "Serum") %>% select(Bile_Acids) %>% unlist(),
          B = venn %>% filter(Type == "Fecal") %>% select(Bile_Acids) %>% unlist(), 
          C = venn %>% filter(Mother_Infant == "Mother") %>% select(Bile_Acids) %>% unlist(),
          D = venn %>% filter(Mother_Infant == "Infant") %>% select(Bile_Acids) %>% unlist())
ggVennDiagram(x, category.names = c("Serum", "Fecal", "Mother", "Infant"), label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
rm(x)

## Plot number of bile acids found in each condition
nb <- venn %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Infant") %>% group_by(Days_Post_Birth, Hypothesis, AMP, Experiment, Sample) %>% count(host_subject_id) #%>% filter(n > 100)
ggplot(nb, aes(x = as.factor(Days_Post_Birth), y = n, fill = AMP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(position=position_dodge(width=0.75),aes(group=AMP, color = AMP)) +
  stat_compare_means(aes(label = after_stat(p.format))) +
  theme_minimal() +
  facet_wrap(~Hypothesis)

## Plot total abundance of BAs in each condition
abt <- inner_join(metadata, TIC, by = "Sample") 
abt <- abt %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite))
abt <- abt %>% pivot_longer(cols =c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
abt <- abt %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Infant") %>% group_by(Days_Post_Birth, Hypothesis, AMP, Experiment, Sample) %>% summarise(n = sum(Expression)*100)
ggplot(abt, aes(x = as.factor(Days_Post_Birth), y = n, fill = AMP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(position=position_dodge(width=0.75),aes(group=AMP, color = AMP)) +
  stat_compare_means(aes(label = after_stat(p.format))) +
  theme_minimal() +
  facet_wrap(~Hypothesis)


## Look at detection frequency for each bile type
total <- venn %>% filter(Type == "Fecal" & Mother_Infant == "Mother") %>% count(Sample)
types <- venn %>% filter(Type == "Fecal" & Mother_Infant == "Mother") %>% group_by(Bile_Type) %>% count(LibraryID) %>% summarise(per = mean(n)/nrow(total)) %>% arrange(Bile_Type)
types$Bile_Type <- factor(types$Bile_Type, levels = c("Penta", "Tetra", "Tri","Di", "Mono", "Non"))
ggplot(types, aes(x = Bile_Type, y = per, fill = Bile_Type)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  coord_flip()

## Look at percent levels of each bile type in blood and feces
types <- venn %>% group_by(Type, Mother_Infant) %>% count(Bile_Type)
ggplot(types, aes(x = Type, y = n, fill = Bile_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Mother_Infant)

## Look at bile acid abundance in fecal samples #####################################################################################
## For this we want to use imputed and batch corrected data 
data <- inner_join(metadata, norm_imp, by = "Sample") 
data <- data %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite))

# Look at bile acid abundance over time
types <- data %>% pivot_longer(cols =c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
types <- inner_join(types, bile_lib, by = c("Bile_Acids" = "Metabolite"), multiple = "any")
types <- types %>% filter(Type == "Fecal") %>% filter(AMP == "FALSE")
types <- types %>% group_by(Sample, Bile_Type, Days_Post_Birth, Mother_Infant) %>% summarise(Total = mean(Expression))
ggplot(types, aes(x=Days_Post_Birth, y=Total, color=Bile_Type)) +
  geom_smooth(method = "loess") +
  theme_minimal() +
  geom_point(alpha = 0.5) +
  facet_wrap(~Mother_Infant, scales = "free_x")

# Look at bile acids with RT matches
RT <- data %>% pivot_longer(cols =c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
RT <- inner_join(RT, bile_lib, by = c("Bile_Acids" = "Metabolite"), multiple = "any")
RT <- RT %>% filter(Type == "Fecal") %>% filter(str_detect(LibraryID, "RTM")) 
RT$Experiment <- str_replace_all(RT$Experiment, ".*FALSE", "Control")
RT$LibraryID <- str_replace_all(RT$LibraryID, ".*RTM", "RTM") 
RT$LibraryID <- str_replace_all(RT$LibraryID, "id.*", "id") 
ggplot(RT, aes(x=Days_Post_Birth, y=Expression, color=as.factor(Experiment))) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  facet_wrap(Mother_Infant~LibraryID, ncol = 5, scales = "free_y")

# Average for each group
average <- data %>% 
  filter(Type == "Fecal") %>% 
  group_by(Days_Post_Birth, Mother_Infant, Type, AMP, Hypothesis) %>% 
  summarise(across(starts_with("X"), ~mean(., na.rm = TRUE)))

# Standard Deviation
sd <- data %>% 
  filter(Type == "Fecal") %>% 
  group_by(Days_Post_Birth, Mother_Infant, Type, AMP, Hypothesis) %>% 
  summarise(across(starts_with("X"), ~sd(., na.rm = TRUE)))

# Fold Change
fold_change <- average %>%
  arrange(Days_Post_Birth, Mother_Infant, Type, Hypothesis)  %>% 
  group_by(Days_Post_Birth, Mother_Infant, Type, Hypothesis) %>% 
  summarise_if(is.numeric, function(x) x[-1]-x[1]) 

#### Significance AMP/Control ##############################################################
### Make sure you're only looking at bile acids that were found in at least 20% of the sample subset
raw <- read.csv("Metabolites_normalized.csv")
raw <- inner_join(metadata, raw, by = "Sample")
raw <- raw %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite))
M <- raw %>% filter(Mother_Infant == "Mother") 
M <- M %>% dplyr::select(where(~sum(. != 0) >= (0.2*nrow(M))))
I <- raw %>% filter(Mother_Infant == "Infant") 
I <- I %>% dplyr::select(where(~sum(. != 0) >= (0.2*nrow(I))))

# Just mom data
pval_M <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Mother")
pval_M <- pval_M %>% filter(Bile_Acids %in% colnames(M))
pval_M <- pval_M %>% group_by(Days_Post_Birth, Mother_Infant, Bile_Acids, Hypothesis) %>% tukey_hsd(Expression ~ as.factor(AMP))
pval_M$p.adj.BH <- p.adjust(pval_M$p.adj, method = "BH")
pval_M$p.adj.char <- ifelse(pval_M$p.adj.BH < 0.001, "***", ifelse(pval_M$p.adj.BH < 0.01, "**", ifelse(pval_M$p.adj.BH < 0.05, "*", "")))
pval_M <- inner_join(pval_M, bile_lib, by = c("Bile_Acids" = "Metabolite"))

pval_I <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Infant")
pval_I <- pval_I %>% filter(Bile_Acids %in% colnames(I))
pval_I <- pval_I %>% group_by(Days_Post_Birth, Mother_Infant, Bile_Acids, Hypothesis) %>% tukey_hsd(Expression ~ as.factor(AMP))
pval_I$p.adj.BH <- p.adjust(pval_I$p.adj, method = "BH")
pval_I$p.adj.char <- ifelse(pval_I$p.adj.BH < 0.001, "***", ifelse(pval_I$p.adj.BH < 0.01, "**", ifelse(pval_I$p.adj.BH < 0.05, "*", "")))
pval_I <- inner_join(pval_I, bile_lib, by = c("Bile_Acids" = "Metabolite"))
rm(raw)

### Line plot of amino acid conjugations 
aa <- pval_I %>% filter(str_detect(LibraryID, "CCM")) %>% filter(str_detect(LibraryID, "Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|Cit|His|Ile|Leu|Lys|Met|Orn|Phe|Pro|Trp|Tyr|Val"))
aa_plot <- data %>% dplyr::select(Sample:Experiment, any_of(aa$Bile_Acids)) %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") 
aa_plot <- aa_plot %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Infant")
aa_plot$Experiment <- str_replace_all(aa_plot$Experiment, ".*FALSE", "Control")

ggplot(aa_plot, aes(x=Days_Post_Birth, y=Expression, color=as.factor(Experiment))) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  scale_x_continuous(limits = c(0,22)) +
  facet_wrap(~Bile_Acids, ncol = 5, scales = "free_y")

##### Heatmap of fold changes - Parent ################################
heat <- fold_change %>%
  filter(Mother_Infant == "Mother") %>%
  ungroup() %>%
  dplyr::select(-Type, -Mother_Infant) %>% 
  dplyr::select(any_of(colnames(M)))
heat$Sample <- paste0(heat$Hypothesis, "_", heat$Days_Post_Birth)
heat$Sample <- factor(heat$Sample, levels = c("1_-3", "1_-2", "1_-1", "1_1", "1_2", "1_7", "1_14", "1_21",
                                              "2_1", "2_2", "2_3", "2_7", "2_14", "2_21", 
                                              "3_1", "3_10", "3_21"))
column_ant <- heat %>% dplyr::select(Sample, Hypothesis) %>% column_to_rownames(var = "Sample")
heat <- heat %>% dplyr::select(-Hypothesis, -Days_Post_Birth) %>% dplyr::select(Sample, everything())
heat <- heat %>%
  gather(Metabolite, value, 2:ncol(heat)) %>%
  spread(Sample, value)
row_ant <- bile_lib %>% filter(Metabolite %in% heat$Metabolite) %>% dplyr::select(Metabolite, Bile_Type) %>% unique()
rownames(row_ant) <- NULL
row_ant <- row_ant %>% column_to_rownames(var = "Metabolite")
heat <- column_to_rownames(heat, var = "Metabolite")

sig_ant <- pval_M %>% dplyr::select(Hypothesis, Days_Post_Birth, Bile_Acids, p.adj.char) %>% pivot_wider(names_from = "Bile_Acids", values_from = "p.adj.char")
sig_ant$Sample <- paste0(sig_ant$Hypothesis, "_", sig_ant$Days_Post_Birth)
sig_ant <- sig_ant %>% dplyr::select(-Hypothesis, -Days_Post_Birth) %>% dplyr::select(Sample, everything())
sig_ant$Sample <- factor(sig_ant$Sample, levels = c("1_-3", "1_-2", "1_-1", "1_1", "1_2", "1_7", "1_14", "1_21",
                                                    "2_1", "2_2", "2_3", "2_7", "2_14", "2_21", 
                                                    "3_1", "3_10", "3_21"))
sig_ant <- sig_ant %>%
  gather(Metabolite, value, 2:ncol(sig_ant)) %>%
  spread(Sample, value)
sig_ant <- column_to_rownames(sig_ant, var = "Metabolite")

breaklist <- seq(-5,5, by = 0.1)
pheatmap(heat, color = colorRampPalette(rev(brewer.pal(7, name = "RdBu")))(length(breaklist)),
         annotation_col = column_ant,
         annotation_row = row_ant,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         cutree_rows = 5, 
         display_numbers = sig_ant, 
         fontsize_number = 6, 
         gaps_col = c(8, 14),
         breaks = breaklist
         
)


##### Heatmap of fold changes - Infant ################################
heat <- fold_change %>%
  filter(Mother_Infant == "Infant") %>%
  ungroup() %>%
  dplyr::select(-Type, -Mother_Infant) %>% 
  dplyr::select(any_of(colnames(I)))
heat$Sample <- paste0(heat$Hypothesis, "_", heat$Days_Post_Birth)
heat$Sample <- factor(heat$Sample, levels = c("1_7", "1_14", "1_21",
                                              "2_7", "2_14", "2_21", 
                                              "3_10", "3_21"))
column_ant <- heat %>% dplyr::select(Sample, Hypothesis) %>% column_to_rownames(var = "Sample")
heat <- heat %>% dplyr::select(-Hypothesis, -Days_Post_Birth) %>% dplyr::select(Sample, everything())
heat <- heat %>%
  gather(Metabolite, value, 2:ncol(heat)) %>%
  spread(Sample, value)
row_ant <- bile_lib %>% filter(Metabolite %in% heat$Metabolite) %>% dplyr::select(Metabolite, Bile_Type) %>% unique()
rownames(row_ant) <- NULL
row_ant <- row_ant %>% column_to_rownames(var = "Metabolite")
heat <- column_to_rownames(heat, var = "Metabolite")

sig_ant <- pval_I %>% filter(Mother_Infant == "Infant") %>% dplyr::select(Hypothesis, Days_Post_Birth, Bile_Acids, p.adj.char) %>% pivot_wider(names_from = "Bile_Acids", values_from = "p.adj.char")
sig_ant$Sample <- paste0(sig_ant$Hypothesis, "_", sig_ant$Days_Post_Birth)
sig_ant <- sig_ant %>% dplyr::select(-Hypothesis, -Days_Post_Birth) %>% dplyr::select(Sample, everything())
sig_ant$Sample <- factor(sig_ant$Sample, levels = c("1_7", "1_14", "1_21",
                                                    "2_7", "2_14", "2_21", 
                                                    "3_10", "3_21"))
sig_ant <- sig_ant %>%
  gather(Metabolite, value, 2:ncol(sig_ant)) %>%
  spread(Sample, value)
sig_ant <- column_to_rownames(sig_ant, var = "Metabolite")

breaklist <- seq(-3,3, by = 0.1)
pheatmap(heat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         annotation_col = column_ant,
         annotation_row = row_ant,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         cutree_rows = 4, 
         display_numbers = sig_ant, 
         fontsize_number = 6, 
         gaps_col = c(3, 6),
         breaks = breaklist
         
)

############## Significance for different days ###########################################
# Mom_data
pdays_M <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Mother_Infant == "Mother")
pdays_M <- pdays_M %>% filter(Bile_Acids %in% colnames(M))
pdays_M <- pdays_M %>% group_by(AMP, Bile_Acids, Hypothesis) %>% tukey_hsd(Expression ~ as.factor(Days_Post_Birth))
pdays_M$p.adj.BH <- p.adjust(pdays_M$p.adj, method = "BH")
pdays_M$p.adj.char <- ifelse(pdays_M$p.adj.BH < 0.001, "***", ifelse(pdays_M$p.adj.BH < 0.01, "**", ifelse(pdays_M$p.adj.BH < 0.05, "*", "")))
pdays_M <- inner_join(pdays_M, bile_lib, by = c("Bile_Acids" = "Metabolite"))

# Infant_data
pdays_I <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Mother_Infant == "Infant")
pdays_I <- pdays_I %>% filter(Bile_Acids %in% colnames(I))
pdays_I <- pdays_I %>% group_by(AMP, Bile_Acids, Hypothesis) %>% tukey_hsd(Expression ~ as.factor(Days_Post_Birth))
pdays_I$p.adj.BH <- p.adjust(pdays_I$p.adj, method = "BH")
pdays_I$p.adj.char <- ifelse(pdays_I$p.adj.BH < 0.001, "***", ifelse(pdays_I$p.adj.BH < 0.01, "**", ifelse(pdays_I$p.adj.BH < 0.05, "*", "")))
pdays_I <- inner_join(pdays_I, bile_lib, by = c("Bile_Acids" = "Metabolite"))

#### Heatmap of average values - Mother ################################
heat <- average %>%
  filter(Mother_Infant == "Mother") %>%
  ungroup() %>%
  dplyr::select(-Type, -Mother_Infant)  %>% 
  dplyr::select(any_of(colnames(M)))
heat$Sample <- paste0(heat$Hypothesis, "_", heat$AMP, "_", heat$Days_Post_Birth)
heat$Sample <- factor(heat$Sample, levels = c("1_FALSE_-3", "1_FALSE_-2", "1_FALSE_-1", "1_FALSE_1", "1_FALSE_2", "1_FALSE_7", "1_FALSE_14", "1_FALSE_21",
                                              "1_TRUE_-3", "1_TRUE_-2", "1_TRUE_-1", "1_TRUE_1", "1_TRUE_2", "1_TRUE_7", "1_TRUE_14", "1_TRUE_21", 
                                              "2_FALSE_1", "2_FALSE_2", "2_FALSE_3", "2_FALSE_7", "2_FALSE_14", "2_FALSE_21",
                                              "2_TRUE_1", "2_TRUE_2", "2_TRUE_3", "2_TRUE_7", "2_TRUE_14", "2_TRUE_21", 
                                              "3_FALSE_1", "3_FALSE_10", "3_FALSE_21",
                                              "3_TRUE_1", "3_TRUE_10", "3_TRUE_21"))

column_ant <- heat %>% dplyr::select(Sample, Days_Post_Birth, AMP, Hypothesis) %>% column_to_rownames(var = "Sample")
heat <- heat %>% dplyr::select(-Days_Post_Birth, -AMP, -Hypothesis) %>% dplyr::select(Sample, everything())
heat <- heat %>%
  gather(Metabolite, value, 2:ncol(heat)) %>%
  spread(Sample, value)
row_ant <- bile_lib %>% filter(Metabolite %in% heat$Metabolite) %>% dplyr::select(Metabolite, Bile_Type) %>% unique()
rownames(row_ant) <- NULL
row_ant <- row_ant %>% column_to_rownames(var = "Metabolite")
heat <- column_to_rownames(heat, var = "Metabolite")

sig_ant <- pdays_M %>% dplyr::select(Hypothesis, AMP, Bile_Acids, group1, group2, p.adj.char) %>% filter((group1 == "-3" & Hypothesis == "1") | (group1 == "1" & Hypothesis != "1"))
sig_ant$Sample <- paste0(sig_ant$Hypothesis, "_", sig_ant$AMP, "_", sig_ant$group2)
sig_ant <- sig_ant %>% dplyr::select(Sample, Bile_Acids, p.adj.char) %>% pivot_wider(names_from = "Bile_Acids", values_from = "p.adj.char")
sig_ant <- sig_ant %>%
  gather(Metabolite, value, 2:ncol(sig_ant)) %>%
  spread(Sample, value)
sig_ant <- column_to_rownames(sig_ant, var = "Metabolite")
sig_ant$"1_FALSE_-3" <- NA
sig_ant$"1_TRUE_-3" <- NA
sig_ant$"2_FALSE_1" <- NA
sig_ant$"2_TRUE_1" <- NA
sig_ant$"3_FALSE_1" <- NA
sig_ant$"3_TRUE_1" <- NA
sig_ant[is.na(sig_ant)] <- ""
sig_ant <- sig_ant %>% dplyr::select("1_FALSE_-3", "1_FALSE_-2", "1_FALSE_-1", "1_FALSE_1", "1_FALSE_2", "1_FALSE_7", "1_FALSE_14", "1_FALSE_21",
                                     "1_TRUE_-3", "1_TRUE_-2", "1_TRUE_-1", "1_TRUE_1", "1_TRUE_2", "1_TRUE_7", "1_TRUE_14", "1_TRUE_21", 
                                     "2_FALSE_1", "2_FALSE_2", "2_FALSE_3", "2_FALSE_7", "2_FALSE_14", "2_FALSE_21",
                                     "2_TRUE_1", "2_TRUE_2", "2_TRUE_3", "2_TRUE_7", "2_TRUE_14", "2_TRUE_21", 
                                     "3_FALSE_1", "3_FALSE_10", "3_FALSE_21",
                                     "3_TRUE_1", "3_TRUE_10", "3_TRUE_21")



breaklist <- seq(-7,7, by = 0.1)
pheatmap(heat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         annotation_col = column_ant,
         annotation_row = row_ant,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         cutree_rows = 4, 
         display_numbers = sig_ant, 
         fontsize_number = 6, 
         gaps_col = c(8,16,22,28,31),
         breaks = breaklist
         
)

#### Heatmap of average values - Infant ################################
heat <- average %>%
  filter(Mother_Infant == "Infant") %>%
  ungroup() %>%
  dplyr::select(-Type, -Mother_Infant) %>% 
  dplyr::select(any_of(colnames(I)))
heat$Sample <- paste0(heat$Hypothesis, "_", heat$AMP, "_", heat$Days_Post_Birth)
heat$Sample <- factor(heat$Sample, levels = c("1_FALSE_7", "1_FALSE_14", "1_FALSE_21",
                                              "1_TRUE_7", "1_TRUE_14", "1_TRUE_21", 
                                              "2_FALSE_7", "2_FALSE_14", "2_FALSE_21",
                                              "2_TRUE_7", "2_TRUE_14", "2_TRUE_21", 
                                              "3_FALSE_10", "3_FALSE_21", "3_TRUE_10", "3_TRUE_21"))


column_ant <- heat %>% dplyr::select(Sample, Days_Post_Birth, AMP, Hypothesis) %>% column_to_rownames(var = "Sample")
heat <- heat %>% dplyr::select(-Days_Post_Birth, -AMP, -Hypothesis) %>% dplyr::select(Sample, everything())
heat <- heat %>%
  gather(Metabolite, value, 2:ncol(heat)) %>%
  spread(Sample, value)
row_ant <- bile_lib %>% filter(Metabolite %in% heat$Metabolite) %>% dplyr::select(Metabolite, Bile_Type) %>% unique()
rownames(row_ant) <- NULL
row_ant <- row_ant %>% column_to_rownames(var = "Metabolite")
heat <- column_to_rownames(heat, var = "Metabolite")

sig_ant <- pdays_I %>% dplyr::select(Hypothesis, AMP, Bile_Acids, group1, group2, p.adj.char) %>% filter(group1 == 7| group1 == 10)
sig_ant$Sample <- paste0(sig_ant$Hypothesis, "_", sig_ant$AMP, "_", sig_ant$group2)
sig_ant <- sig_ant %>% dplyr::select(Sample, Bile_Acids, p.adj.char) %>% pivot_wider(names_from = "Bile_Acids", values_from = "p.adj.char")
sig_ant <- sig_ant %>%
  gather(Metabolite, value, 2:ncol(sig_ant)) %>%
  spread(Sample, value)
sig_ant <- column_to_rownames(sig_ant, var = "Metabolite")
sig_ant$"1_FALSE_7" <- NA
sig_ant$"1_TRUE_7" <- NA
sig_ant$"2_FALSE_7" <- NA
sig_ant$"2_TRUE_7" <- NA
sig_ant$"3_FALSE_10" <- NA
sig_ant$"3_TRUE_10" <- NA
sig_ant[is.na(sig_ant)] <- ""
sig_ant <- sig_ant %>% dplyr::select("1_FALSE_7", "1_FALSE_14", "1_FALSE_21",
                                     "1_TRUE_7", "1_TRUE_14", "1_TRUE_21", 
                                     "2_FALSE_7", "2_FALSE_14", "2_FALSE_21",
                                     "2_TRUE_7", "2_TRUE_14", "2_TRUE_21", 
                                     "3_FALSE_10", "3_FALSE_21", "3_TRUE_10", "3_TRUE_21")

breaklist <- seq(-6,6, by = 0.1)
pheatmap(heat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         annotation_col = column_ant,
         annotation_row = row_ant,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         cutree_rows = 4, 
         display_numbers = sig_ant, 
         fontsize_number = 6, 
         gaps_col = c(3,6,9, 12, 14),
         breaks = breaklist
         
)

####### Significance comparing hypotheses #################################################
phyp_M <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Mother_Infant == "Mother") %>% filter(Days_Post_Birth == 21)
phyp_M <- phyp_M %>% filter(Bile_Acids %in% colnames(M))
phyp_M <- phyp_M %>% group_by(Bile_Acids, AMP) %>% tukey_hsd(Expression ~ as.factor(Hypothesis))
phyp_M$p.adj.BH <- p.adjust(phyp_M$p.adj, method = "BH")
phyp_M$p.adj.char <- ifelse(phyp_M$p.adj.BH < 0.001, "***", ifelse(phyp_M$p.adj.BH < 0.01, "**", ifelse(phyp_M$p.adj.BH < 0.05, "*", "")))

phyp_I <- data %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") %>% filter(Mother_Infant == "Infant") %>% filter(Days_Post_Birth == 21)
phyp_I <- phyp_I %>% filter(Bile_Acids %in% colnames(I))
phyp_I <- phyp_I %>% group_by(Bile_Acids, AMP) %>% tukey_hsd(Expression ~ as.factor(Hypothesis))
phyp_I$p.adj.BH <- p.adjust(phyp_I$p.adj, method = "BH")
phyp_I$p.adj.char <- ifelse(phyp_I$p.adj.BH < 0.001, "***", ifelse(phyp_I$p.adj.BH < 0.01, "**", ifelse(phyp_I$p.adj.BH < 0.05, "*", "")))
phyp_I <- inner_join(bile_lib, phyp_I, by = c("Metabolite" = "Bile_Acids"))

## Line plot
b <- "X2546_i"

plot <- data %>% dplyr::select(Sample:Experiment, !!sym(b)) %>% filter(Mother_Infant == "Infant")

ggplot(plot, aes(x=Days_Post_Birth, y=!!sym(b), color=as.factor(Hypothesis), shape=as.factor(Hypothesis))) +
  geom_smooth(method = "loess") +
  geom_point() +
  theme_minimal() +
  facet_wrap(~AMP, scales = "free_x") 

############ Volcano Plot ####################################
volcano <- inner_join(metadata, norm_imp, by = "Sample") 
volcano <- volcano %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite))
volcano <- volcano %>% filter(Type == "Fecal") %>% filter(Mother_Infant == "Mother") %>% 
                       filter((Days_Post_Birth == 21 & Hypothesis == 2)|(Days_Post_Birth == 21 & Hypothesis == 2))
## Only include compounds that were found in 50% of samples
test <- norm %>% dplyr::select(Sample, any_of(colnames(volcano))) %>% filter(Sample %in% volcano$Sample)
test[test == 0] <- NA
test <- test %>% dplyr::select(where(~sum(is.na(.)) >= (0.5*nrow(test))))
volcano <- volcano %>% dplyr::select(!any_of(colnames(test)))

sig <- volcano %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") 
sig <- sig %>% group_by(Bile_Acids) %>% tukey_hsd(Expression ~ as.factor(Hypothesis))
sig$p.adj.BH <- p.adjust(sig$p.adj, method = "BH")
sig$p.adj.char <- ifelse(sig$p.adj.BH < 0.001, "***", ifelse(sig$p.adj.BH < 0.01, "**", ifelse(sig$p.adj.BH < 0.05, "*", "")))

volcano <- volcano %>% group_by(Hypothesis) %>% summarise(across(starts_with("X"), ~mean(., na.rm = TRUE)))
volcano <- volcano %>% summarise_if(is.numeric, function(x) x[-1]-x[1]) 
volcano <- volcano %>% pivot_longer(starts_with("X"), values_to = "FC", names_to = "Bile_Acids")

volcano <- left_join(volcano, sig, by = "Bile_Acids")
volcano <- inner_join(volcano, bile_lib, by = c("Bile_Acids" = "Metabolite"), multiple = "all")
volcano$ID <- str_replace_all(volcano$LibraryID, c("Candidate " = "", ";.*" = "", " bile acid..delta mass" = "; delta"))
volcano$ID <- ifelse(volcano$p.adj.BH > 0.05,"", volcano$ID)

volcano$dfex <- "NO"
volcano$dfex[volcano$FC > 0.5 & volcano$p.adj.BH < 0.05] <- "UP"
volcano$dfex[volcano$FC < -0.5 & volcano$p.adj.BH < 0.05] <- "DOWN"

volcano[is.na(volcano)] <- ""

# Volcano Plot
ggplot(volcano, aes(x=FC, y=-log10(p.adj.BH), col=dfex, label=ID)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept=c(-0.5, 0.5)) +
  geom_hline(yintercept=-log10(0.05))

## PCA and PLS-DA ###########################################################################
pd <- norm_imp %>% dplyr::select(Sample, any_of(bile_lib$Metabolite))
pdm <- metadata %>% filter(Mother_Infant == "Infant") %>% filter(Type == "Fecal") %>% filter(Days_Post_Birth == 21) %>% filter(Hypothesis == 1)
pd <- pd %>% filter(Sample %in% pdm$Sample)
pdm <- pdm %>% filter(Sample %in% pd$Sample)
pd <- pd %>% column_to_rownames("Sample")

PCA <- mixOmics::pca(pd, ncomp = 2, scale = FALSE, center = FALSE)
PCA_scores <- data.frame(PCA$variates$X, pdm)
plotIndiv(PCA, ellipse = TRUE, group = pdm$Experiment, centroid = TRUE, legend = TRUE, ind.names = FALSE, pch = 20)
# Permanova
dist <- vegdist(pd, method = "euclidean")
permanova <- adonis2(dist ~ PCA_scores$AMP, PCA_scores, na.action = na.omit)

PLSDA <- plsda(pd, as.factor(pdm$Experiment), scale = FALSE)
plotIndiv(PLSDA, ellipse = TRUE, group = pdm$Experiment, centroid = TRUE, legend = TRUE, ind.names = FALSE, pch = 20)

perf_PLSDA <- perf(PLSDA, validation = "Mfold", auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA, legend = TRUE)

#plotVar(PLSDA)

VIPs <- as.data.frame(mixOmics::vip(PLSDA))
VIPs <- VIPs %>% rownames_to_column("Metabolite") 
VIPs <- inner_join(VIPs, bile_lib, by = "Metabolite")
VIPs$ID <- str_replace_all(VIPs$LibraryID, c("Candidate " = "", ";.*" = "", " bile acid..delta mass" = "; delta"))
VIPs$ID <- paste0(VIPs$Metabolite, "; ", VIPs$ID)
names <- VIPs$ID %>% as.factor()
plotLoadings(PLSDA, method = 'median', contrib = 'max', ndisplay = 25, name.var = names)
write_csv(VIPs, "VIPs_3.csv")

## Look at weights ##########################################################################
data_w <- inner_join(metadata, norm_imp, by = "Sample") 
data_w <- data_w %>% dplyr::select(Sample:Experiment, any_of(bile_lib$Metabolite)) %>% filter(Mother_Infant == "Infant") %>% filter(Days_Post_Birth > 19) %>% filter(Type == "Fecal")
weights <- read.csv("Weight_Gain.csv")
# Have to get rid of some weights where sex wasn't collected
weights <- weights %>% filter(!is.na(Sex))

# Bloxplot of weights
ggplot(weights, aes(x = Sex, y = Weight, color = AMP)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(label = "p.format") +
  geom_point(position = position_dodge(width=0.75)) +
  facet_wrap(~Hypothesis)

## Spearman correlation of bile acids and weight
weights <- weights %>% group_by(Hypothesis, Batch, AMP, Cage, Sex) %>% summarise(Mean = mean(Weight))
weights$Cage <- as.character(weights$Cage)
data_w <- inner_join(weights, data_w, by = c("Hypothesis", "Batch", "AMP", "Cage")) 
data_w <- data_w %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression") 

corr <- data_w %>% group_by(Sex, Bile_Acids) %>% rstatix::cor_test(vars=c("Mean"), vars2=c("Expression"), method="spearman")
corr$p.adj.BH <- p.adjust(corr$p, method = "BH")
corr$p.adj.char <- ifelse(corr$p.adj.BH < 0.001, "***", ifelse(corr$p.adj.BH < 0.01, "**", ifelse(corr$p.adj.BH < 0.05, "*", "")))
corr <- left_join(corr, bile_lib, by = c("Bile_Acids" = "Metabolite"))

# Graph correlation for specific metabolite
i <- "X27148"

bile_plot <- data_w %>% filter(Bile_Acids == i)
bile_plot$Hypothesis <- as.factor(bile_plot$Hypothesis)
ggscatter(bile_plot, x = "Mean", y = "Expression", shape = "Hypothesis", add = "reg.line",  conf.int = TRUE, 
          cor.method = "spearman") + stat_cor() + facet_wrap(~Sex)

######### Look at bile acids over time using SantaR ##############################################
snt <- read.csv("Metabolites_normalized.csv")
snt[snt == 0] <- NA

snt <- snt %>% dplyr::select(Sample, any_of(bile_lib$Metabolite))
## Choose which variables you want to look at 
meta <- metadata %>% filter(AMP == TRUE) %>% filter(Mother_Infant == "Mother") %>% filter(Hypothesis < 3) %>% filter(Type == "Fecal")
test <- meta %>% count(host_subject_id) %>% filter(n < 4)
meta <- meta %>% filter(!host_subject_id %in% test$host_subject_id)
meta$Corrected_Time <- ifelse(meta$Experiment == "1_TRUE", meta$Days_Post_Birth+4, meta$Days_Post_Birth)

snt <- snt %>% filter(Sample %in% meta$Sample)
meta <- meta %>% filter(Sample %in% snt$Sample)
snt <- snt %>% column_to_rownames(var = "Sample")
meta <- meta %>% column_to_rownames(var = "Sample")

## Only include compounds that were found in 50% of samples
snt <- snt %>% dplyr::select(where(~sum(is.na(.)) <= (0.5*nrow(snt))))

res_AMP <- santaR_auto_fit(inputData=snt, ind=meta$host_subject_id, time=meta$Corrected_Time, group=meta$Hypothesis, df=5, ncores=6, CBand=TRUE, pval.dist=TRUE)

## Generate a summary
pval_res <- santaR_auto_summary(SANTAObjList=res_AMP, targetFolder=NA)
pval_res$pval.summary
pval_all <- pval_res$pval.all
pval_all <- pval_all %>% rownames_to_column("Metabolite")
pval_all <- inner_join(pval_all, bile_lib, by = "Metabolite")

santaR_plot(res_AMP$X37451, showIndPoint=TRUE, showIndCurve=TRUE, showGroupMeanCurve=TRUE, showConfBand=TRUE)

i <- "X27148"
plot <- inner_join(metadata, norm_imp, by = "Sample") 
plot <- plot %>% filter(AMP == TRUE) %>% filter(Mother_Infant == "Mother") %>% filter(Hypothesis < 3) %>% filter(Type == "Fecal")

ggplot(plot, aes(x=Days_Post_Birth, y=!!sym(i), color=as.factor(Experiment))) +
  geom_smooth(method = "loess") +
  geom_point(alpha = 0.5) +
  theme_minimal()

######## Run dimensionality reduction over time (TEMPTED) ####################################
tmpt <- inner_join(metadata, norm_imp, by = "Sample")  %>% filter(Mother_Infant == "Infant") %>% filter(Type == "Fecal")
test <- tmpt %>% count(host_subject_id) %>% filter(n <= 1)
tmpt <- tmpt %>% filter(!host_subject_id %in% test$host_subject_id)

meta <- tmpt %>% dplyr::select(Sample:Experiment)
tmpt <- tmpt %>% dplyr::select(Sample, starts_with("X"))
#tmpt <- tmpt %>% dplyr::select(Sample, any_of(bile_lib$Metabolite))
meta <- meta %>% column_to_rownames("Sample")
tmpt <- tmpt %>% column_to_rownames("Sample")

## Check that rownames match between meta and tmpt
table(rownames(tmpt) == rownames(meta))

## Run Tempted on preprocessed data
res <- tempted_all(tmpt,
                   meta$Days_Post_Birth,
                   meta$host_subject_id,
                   threshold=1,
                   transform="none",
                   r=2,
                   smooth=1e-5,
                   do_ratio = FALSE,
                   pct_aggregate=0.075)

## Plot PCoA
## You can choose which variable to plot
b <- "Experiment"

uni <- meta %>% 
  dplyr::select(host_subject_id, all_of(b), Hypothesis, Batch) %>% 
  unique()
rownames(uni) <- NULL
uni <- uni %>% column_to_rownames("host_subject_id")
uni <- cbind(res$A_hat[rownames(uni),], uni)

## Plot samples
ggplot(data=uni, aes(x=uni[,1], y=uni[,2], color=!!sym(b))) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='subject loading') +
  stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = !!sym(b))) +
  facet_wrap(~Hypothesis)

## Plot principal components over time
plot_time_loading(res, r=2) + 
  geom_line(linewidth=1) + 
  labs(title='temporal loadings', x='days')

## Plot features
features <- as.data.frame(res$B_hat)
features <- rownames_to_column(features, var = "Metabolite") %>% rename(PC1 = 'Component 1', PC2 = 'Component 2')
features <- left_join(features, bile_lib, by = "Metabolite")
ggplot(features, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = Bile_Type)) + 
  labs(x='Component 1', y='Component 2', title='feature loading') 

## Plot group trajectories
group <- meta %>% 
  dplyr::select(host_subject_id, all_of(b)) %>% 
  unique()
plot_metafeature(res$metafeature_aggregate, group) + xlab("Days of Life")

## Find VIPs
VIP <- as.data.frame(res$toppct_aggregate)
VIP <- cbind(features, VIP)
all <- read.csv("library_matches_all.csv")
all <- all %>% filter(!is.na(LibraryID)) %>% dplyr::select(-Network_Number)
all$Metabolite <- paste0("X", all$Metabolite)
VIP <- left_join(VIP, all, by = "Metabolite")
VIP$LibraryID.x <- ifelse(is.na(VIP$LibraryID.x), VIP$LibraryID.y, VIP$LibraryID.x)
VIP$Precursor_Mass.x <- ifelse(is.na(VIP$Precursor_Mass.x), VIP$Precursor_Mass.y, VIP$Precursor_Mass.x)
VIP$RTConsensus.x <- ifelse(is.na(VIP$RTConsensus.x), VIP$RTConsensus.y, VIP$RTConsensus.x)
VIP <- VIP %>% dplyr::select(-ends_with("y")) %>% rename_with(~str_replace_all(., ".x", ""))
rownames(VIP) <- NULL
write_csv(VIP, "VIPs.csv")

## Plot annotated important features
VIPf <- VIP %>% filter(`Component 1` == "TRUE") %>% filter(!is.na(LibraryID))
VIPf$Bile <- ifelse(str_detect(VIPf$LibraryID, "bile"), "Bile", "Other")
ggplot(VIPf, aes(x = reorder(Metabolite, PC1), y = PC1, fill = Bile)) +
  geom_bar(stat = "identity") +
  coord_flip()

## Plot all VIPs for PC1 or PC2
VIP$C1 <- ifelse(VIP$PC1 < 0, "Neg", "Pos")
VIP$C2 <- ifelse(VIP$PC2 < 0, "Neg", "Pos")
VIP1 <- VIP %>% filter(`Component 1` == "TRUE")
VIP1 <- tmpt %>% dplyr::select(any_of(VIP1$Metabolite)) %>% rownames_to_column("Sample")
VIP1 <- inner_join(metadata, VIP1, by = "Sample")
VIP1 <- VIP1 %>% pivot_longer(cols = c(starts_with("X")), names_to = "Bile_Acids", values_to = "Expression")
VIP1 <- inner_join(VIP1, VIP, by = c("Bile_Acids" = "Metabolite"), multiple = "any")
VIP1$LibraryID <- str_remove_all(VIP1$LibraryID, " OR .*")
VIP1$LibraryID <- str_remove_all(VIP1$LibraryID, "Candidate ")
VIP1$LibraryID <- str_remove_all(VIP1$LibraryID, "sirius.*")
VIP1$LibraryID <- str_remove_all(VIP1$LibraryID, "buddy.*")

ggplot(VIP1, aes(x=Days_Post_Birth, y=Expression, color=as.factor(Bile_Acids))) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(C1~Experiment)

## Plot all samples
feat <- res$metafeature_aggregate
colnames(feat)[2] <- 'host_subject_id'
feat <- reshape(feat, 
                idvar=c("host_subject_id","timepoint") , 
                v.names=c("value"), timevar="PC",
                direction="wide")
colnames(feat) <- sub(".*value[.]Component ", "Component_",  colnames(feat))
feat <- left_join(feat, group, by = "host_subject_id", multiple = "all")


ggplot(data=feat, aes(x=Component_1, y=Component_2)) +
  geom_point(aes(color = timepoint)) + scale_color_distiller(palette='RdBu') + 
  labs(x='Component 1', y='Component 2', color='Day')

ggplot(data=feat, aes(x=Component_1, y=Component_2, color=!!sym(b))) +
  geom_point() + 
  labs(x='Component 1', y='Component 2')
