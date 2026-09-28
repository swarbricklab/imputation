##### Read in Libraries #####
library(tidyverse)
library(ggpubr)
library(cluster)
library(RColorBrewer)
library(dplyr)

# Log output and messages
log_file_conn <- file(snakemake@log[[1]], open = "wt")
sink(log_file_conn)
sink(log_file_conn, type = "message")

message("Getting output paths from rule")
sexcheck_path <- snakemake@output[['sexcheck']]
anc_check_path <- snakemake@output[['anc_check']]
plot_path <- snakemake@output[['plot']]
outdir <- dirname(sexcheck_path)

message("Loading data scores")
data_score <- read_tsv(
    snakemake@input[['projected_scores']], 
    col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "Provided_Ancestry" = "c")
  ) %>%
  rename(FID = `#FID`)
print(head(data_score))

message("Loading reference scores")
onekg_score <- read_tsv(
    snakemake@input[['projected_1000g_scores']], 
    col_types = cols(.default = "d", "#IID" = "c", "SuperPop" = "c", "Population" = "c")
  ) %>%
  rename(IID = `#IID`)
print(head(onekg_score))

message("Loading reference ancestry")
onekg_anc <- read_tsv(
    snakemake@input[['fam_1000g']], 
    col_types = cols(.default = "c", "SEX" = "d")
  ) %>%
  rename(IID = `#IID`)
print(head(onekg_anc))

message("Loading data ancestry")
data_anc <- read_tsv(
    snakemake@input[['psam']], 
    col_types = cols(.default = "c", "SEX" = "d")
  ) %>%
  rename(FID = `#FID`) 
print(head(data_anc))

message("Loading sex check information")
sex_check <- read_tsv(
    snakemake@input[['sexcheck']], 
    col_types = cols(.default = "c", "PEDSEX" = "d", "SNPSEX" = "d", "F" = "d")
  ) 
print(head(sex_check))

# Select columns IID and those that start with "PC" using select() and everything()
onekg_score_temp <- onekg_score %>%
  select(IID, starts_with("PC")) %>%
  rename_with(~ str_remove(., "_AVG$")) %>%
  mutate(FID = NA) %>%
  select(FID, IID, everything())

message("1000g score:")
print(head(onekg_score_temp))

data_score_temp <- data_score %>%
  select(any_of(c("FID")), IID, starts_with("PC")) %>%
  rename_with(~ str_remove(., "_AVG$")) %>%
  select(FID, IID, everything())

message("Data score:")
print(head(data_score_temp))

colnames(onekg_score_temp)
colnames(data_score_temp)

scores <- rbind(onekg_score_temp, data_score_temp)
print(head(as.data.frame(scores)))

scores <- left_join(scores, onekg_anc %>% select(IID, SuperPop))
scores <- left_join(scores, data_anc %>% select(IID, Provided_Ancestry))
print(head(scores))

##### Calculate Medoids and Assign Clusters #####
# Cluster on the 10 principal components only. NB: positional `scores[,2:11]`
# wrongly included the IID (character) column and dropped PC10; R 3.6 tolerated
# the stray character column but R >= 4 clusters on noise, collapsing every
# sample to the largest (EUR) cluster. Select the PC columns explicitly.
pc_cols <- paste0("PC", 1:10)
pam_res <- pam(scores[, pc_cols], 6)

##### Assign Ancestries to Individuals #####
scores$Cluster <- factor(pam_res$clustering)
print(head(as.data.frame(scores)))

conversion_table <- table(scores$SuperPop, scores$Cluster)
conversion_key <- data.frame(Cluster = colnames(conversion_table), Assignment = NA)

for (clust in conversion_key$Cluster){
  conversion_key$Assignment[which(conversion_key$Cluster == clust)] <- rownames(conversion_table)[which.max(conversion_table[,clust])]
}
print(conversion_key)

scores <- left_join(scores, conversion_key)

scores$combined_assignment <- ifelse(is.na(scores$SuperPop), scores$Assignment, scores$SuperPop)
scores$combined_assignment <- ifelse(is.na(scores$SuperPop), scores$Assignment, scores$SuperPop)
scores$Changed <- ifelse(is.na(scores$Provided_Ancestry), "Matched", ifelse(scores$Provided_Ancestry == scores$Assignment, "Matched", paste0("Unmatched-",scores$Provided_Ancestry,"->",scores$Assignment)))

##### Plot Results #####
df4plots <- rbind(data.frame(scores[which(is.na(scores$Provided_Ancestry)),], Plot = "1000G Reference"), data.frame(scores[which(!is.na(scores$Provided_Ancestry)),], Plot = "Projected Data Assignments"), data.frame(scores[which(!is.na(scores$Provided_Ancestry)),], Plot = "Projected Data Assignments vs Original Assignments"))
df4plots$Final_Assignment <- ifelse(df4plots$Plot == "Projected Data Assignments", df4plots$Assignment, ifelse(df4plots$Plot == "Projected Data Assignments vs Original Assignments", df4plots$Changed, df4plots$combined_assignment))
df4plots <- arrange(df4plots, Final_Assignment)
df4plots$Final_Assignment <- factor(df4plots$Final_Assignment, levels = c(unique(df4plots$Assignment), unique(df4plots$Changed)))

##### Set up population colors #####
matching_colors <- c("gray88",colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(unique(scores$Changed))-1))
names(matching_colors) <- sort(unique(scores$Changed))

pop_colors <- brewer.pal(length(unique(scores$Assignment)), "Dark2")
names(pop_colors) <- unique(scores$Assignment)

colors <- c(matching_colors, pop_colors)

plot_PCs_medoids <- ggplot(df4plots, aes(PC1, PC2, color = Final_Assignment)) +
  geom_point() +
  theme_bw() +
  facet_wrap(vars(Plot)) +
  scale_color_manual(values = colors)

ggsave(plot_PCs_medoids, filename = plot_path, height = 5, width = 12)

##### Subset the Mismatching Ancestry and Sex data #####

sex_check %>%
  filter(STATUS == "PROBLEM") %>%
  mutate(`UPDATE/REMOVE/KEEP` = "UPDATE") %>%
  write_tsv(sexcheck_path, na = "")

anc_mismatch <- df4plots[which(df4plots$Plot == "Projected Data Assignments vs Original Assignments" & df4plots$Final_Assignment != "Matched"),]
anc_temp <- df4plots[,c("FID","IID", "Assignment")]
colnames(anc_temp) <- c("FID","IID", "PCA_Assignment")

anc_mismatch <- left_join(anc_mismatch[,c("FID","IID")], data_anc)
anc_mismatch <- left_join(anc_mismatch, anc_temp)
anc_mismatch <- unique(anc_mismatch)

anc_mismatch <- anc_mismatch %>%
  select("FID", "IID", "Provided_Ancestry", "PCA_Assignment") %>%
  mutate(`UPDATE/REMOVE/KEEP` = "UPDATE") %>%
  rename(`#FID` = FID) %>%
  write_tsv(anc_check_path, na = "")

sink()