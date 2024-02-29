##### Project: Osmia Climate Change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of ITS data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(stringr) # Version 1.5.1  
  library(ggplot2) # Version 3.4.3
  library(phyloseq) # Version 1.44.0
  library(plotrix) # Version 3.8-4
  library(vegan) # Version 2.6-4
  library(magrittr) # Version 2.0.3
  library(tidyverse) # Version 1.2.0
  library(decontam) # Version 1.20.0
  library(microbiome) # Version 1.22.0
  library(car) # Version 3.1-2
  library(nlme) # Version 3.1-164
  library(emmeans) # Version 1.10.0
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(DESeq2) # Version 1.40.2

# Italicize title
  fung.title <- expression(paste("(", italic("b"), ")"))

# Set color scheme
  climate.colors <- c("CS" = "#64B5F6",
                      "CN" = "#1565C0",
                      "AS" = "#9E9E9E",
                      "AN" = "#616161",
                      "WS" = "#E57373",
                      "WN" = "#C62828")

# Set labels
  climate.labs <- c("CS" = "Cool: Sterile",
                    "CN" = "Cool: Natural",
                    "AS" = "Ambient: Sterile",
                    "AN" = "Ambient: Natural",
                    "WS" = "Warm: Sterile",
                    "WN" = "Warm: Natural") 

# Import data
  seqtab.nochim1 <- readRDS("OsmiaCC_seqsITS_run1.rds")
  taxa1 <- readRDS("OsmiaCC_taxaITS_run1.rds")
  metaITS.CC.run1 <- read.csv("OsmiaCC_master - ITS_run1.csv")
  
  seqtab.nochim2 <- readRDS("OsmiaCC_seqsITS_run2.rds")
  taxa2 <- readRDS("OsmiaCC_taxaITS_run2.rds")
  metaITS.CC.run2 <- read.csv("OsmiaCC_master - ITS_run2.csv")

## Create phyloseq objects for each ITS run ----

# Re-create your df
  samples.out1 <- rownames(seqtab.nochim1)
  samples.out1 <- stringr::str_sort(samples.out1, numeric = TRUE)
  samples1 <- data.frame(metaITS.CC.run1)
  extractionID <- samples1$extractionID
  sample_type <- samples1$sample_type
  sampleID <- samples1$sampleID
  temp_treat <- samples1$temp_treat
  micro_treat <- samples1$micro_treat
  combo_treat <- samples1$combo_treat
  sample_or_control <- samples1$sample_or_control
  sex <- samples1$sex
  graft_stage <- samples1$graft_stage
  DNA_conc <- samples1$DNA_conc
  sampleinfo1 <- data.frame(extractionID = extractionID, 
                            sample_type = sample_type,
                            sampleID = sampleID,  
                            temp_treat = temp_treat, 
                            micro_treat = micro_treat, 
                            combo_treat = combo_treat, 
                            sample_or_control = sample_or_control,
                            sex = sex,
                            graft_stage = graft_stage,
                            DNA_conc = DNA_conc)
  rownames(sampleinfo1) <- samples.out1

# Format your data to work with phyloseq
  ps8 <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows = FALSE), sample_data(sampleinfo1), tax_table(taxa1))
  ps8
  
# Re-create your df
  samples.out2 <- rownames(seqtab.nochim2)
  samples.out2 <- stringr::str_sort(samples.out2, numeric = TRUE)
  samples2 <- data.frame(metaITS.CC.run2)
  extractionID <- samples2$extractionID
  sample_type <- samples2$sample_type
  sampleID <- samples2$sampleID
  temp_treat <- samples2$temp_treat
  micro_treat <- samples2$micro_treat
  combo_treat <- samples2$combo_treat
  sample_or_control <- samples2$sample_or_control
  sex <- samples2$sex
  graft_stage <- samples2$graft_stage
  DNA_conc <- samples2$DNA_conc
  sampleinfo2 <- data.frame(extractionID = extractionID,
                            sample_type = sample_type,
                            sampleID = sampleID,
                            temp_treat = temp_treat,
                            micro_treat = micro_treat,
                            combo_treat = combo_treat,
                            sample_or_control = sample_or_control,
                            sex = sex,
                            graft_stage = graft_stage,
                            DNA_conc = DNA_conc)
  rownames(sampleinfo2) <- samples.out2
  
# Format your data to work with phyloseq
  ps9 <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows = FALSE), sample_data(sampleinfo2), tax_table(taxa2))
  ps9
  
# Merge phyloseq objects
  ps10 <- merge_phyloseq(ps8, ps9)
  ps10
  
# Display total number of reads, mean, and se in phyloseq obj before processing
  sum(sample_sums(ps10))
  mean(sample_sums(ps10))
  print(plotrix::std.error(sample_sums(ps10)))
  
## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Create df with LibrarySize for each sample
  df <- as.data.frame(sample_data(ps10))
  df$LibrarySize <- sample_sums(ps10)
  df <- df[order(df$LibrarySize), ]
  df$Index <- seq(nrow(df))
  
# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_or_control)) + 
    geom_point()
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps10)$is.neg <- sample_data(ps10)$sample_or_control == "control"
  contamdf.prev <- decontam::isContaminant(ps10, method = "prevalence", neg = "is.neg", threshold = 0.1)
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls with a higher threshold
  contamdf.prev05 <- decontam::isContaminant(ps10, method = "prevalence", neg = "is.neg", threshold = 0.5)
  table(contamdf.prev05$contaminant)
  head(which(contamdf.prev05$contaminant))
  
# Make phyloseq object of presence-absence in negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps10)$sample_or_control == "control", ps10)
  
# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))
  
# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps10)$sample_or_control == "sample", ps10)
  
# Calculate taxa abundance in samples from sample counts
  ps.pos.presence <- phyloseq::transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))
  
# Make data.frame of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contam.prev05 = contamdf.prev05$contaminant)
  
# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.prev05)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)") 
  
# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.prev05$contaminant, ps10)
  ps.noncontam
  
# Remove control samples used for identifying contaminants
  ps.sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps.sub
  
# Remove samples without any reads
  ps11 <- phyloseq::prune_samples(sample_sums(ps.sub) != 0, ps.sub)
  ps11
  
# Remove patterns in tax_table   
  tax_table(ps11)[, colnames(tax_table(ps11))] <- gsub(tax_table(ps11)[, colnames(tax_table(ps11))], pattern = "[a-z]__", replacement = "")  
  
# Transform counts to relative abundances
  ps11.relabund <- phyloseq::transform_sample_counts(ps11, function(x) x/sum(x))
  
# Save taxonomy, raw reads, and relative abundance data
  write.csv(tax_table(ps11), "OsmiaCC_ITStaxa_all.csv")
  write.csv(otu_table(ps11), "OsmiaCC_ITSotu_all.csv")
  write.csv(otu_table(ps11.relabund), "OsmiaCC_ITS_relabund.csv")
  
# Subset provisions collected before and after homogenization
  ps12 <- phyloseq::subset_samples(ps11, sample_type == "initial provision")
  ps12
  
# Provisions without bees
  
# Subset daata  
  ps13 <- phyloseq::subset_samples(ps11, sample_type == "provision w/o bee")
  ps13
  
# Display total number of reads, mean, and se in phyloseq obj after processing
  sum(sample_sums(ps13))
  mean(sample_sums(ps13))
  print(plotrix::std.error(sample_sums(ps13)))
  
# Calculate the reads per sample
  reads.sample.NoBee <- microbiome::readcount(ps13)
  head(reads.sample.NoBee)
  
# Add reads per sample to meta data
  sample_data(ps13)$reads.sample.NoBee <- reads.sample.NoBee
  
# Save sample metadata
  meta.NoBee <- sample_data(ps13)
  
# How many samples for each developmental stage?
  meta.NoBee %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.NoBee),
              se = sd(reads.sample.NoBee)/sqrt(N),
              max = max(reads.sample.NoBee),
              min = min(reads.sample.NoBee))
  
# Add Seq to each taxa name
  taxa_names(ps13) <- paste0("Seq", seq(ntaxa(ps13)))
  
# Create a df containing the number of reads per OTU
  read.sums.df.NoBee <- data.frame(nreads = sort(taxa_sums(ps13), TRUE), 
                           sorted = 1:ntaxa(ps13),
                           type = "OTUs")
  
# Add a column containing the number of reads per sample
  read.sums.df.NoBee <- rbind(read.sums.df.NoBee, data.frame(nreads = sort(sample_sums(ps13), TRUE), 
                                                             sorted = 1:nsamples(ps13),
                                                             type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(read.sums.df.NoBee, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") +
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
# Provisions with bees
  
# Subset daata  
  ps14 <- phyloseq::subset_samples(ps11, sample_type == "final provision")
  ps14
  
# Display total number of reads, mean, and se in phyloseq obj after processing
  sum(sample_sums(ps14))
  mean(sample_sums(ps14))
  print(plotrix::std.error(sample_sums(ps14)))
  
# Calculate the reads per sample
  reads.sample.bee <- microbiome::readcount(ps14)
  head(reads.sample.bee)
  
# Add reads per sample to meta data
  sample_data(ps14)$reads.sample.bee <- reads.sample.bee
  
# Save sample metadata
  meta.bee <- sample_data(ps14)
  
# How many samples for each developmental stage?
  meta.bee %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.bee),
              se = sd(reads.sample.bee)/sqrt(N),
              max = max(reads.sample.bee),
              min = min(reads.sample.bee))

# Subset data by sex
  meta.bee.M <- meta.bee[meta.bee$sex == "M", ]
  meta.bee.F <- meta.bee[meta.bee$sex == "F", ]
    
# Read stats for provisions with bees - males
  meta.bee.M %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.bee),
              se = sd(reads.sample.bee)/sqrt(N),
              max = max(reads.sample.bee),
              min = min(reads.sample.bee))
  
# Read stats for provisions with bees - males
  meta.bee.F %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.bee),
              se = sd(reads.sample.bee)/sqrt(N),
              max = max(reads.sample.bee),
              min = min(reads.sample.bee))
  
# Add Seq to each taxa name
  taxa_names(ps14) <- paste0("Seq", seq(ntaxa(ps14)))
  
# Create a df containing the number of reads per OTU
  read.sums.df.bee <- data.frame(nreads = sort(taxa_sums(ps14), TRUE), 
                                 sorted = 1:ntaxa(ps14),
                                 type = "OTUs")
  
# Add a column containing the number of reads per sample
  read.sums.df.bee <- rbind(read.sums.df.bee, data.frame(nreads = sort(sample_sums(ps14), TRUE), 
                                                         sorted = 1:nsamples(ps14),
                                                         type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(read.sums.df.bee, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") +
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
# Provisions with and without bees
  
# Subset data
  ps15 <- phyloseq::subset_samples(ps11, sample_type != "initial provision")
  ps15
  
## Richness and alpha diversity ----
  
# Provisions without bees  
  
# Estimate richness and alpha diversity
  fung.rich.NoBee <- phyloseq::estimate_richness(ps13, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  fung.rich.NoBee$sampleID <- sample_data(ps13)$sampleID
  fung.rich.NoBee$sample_type <- sample_data(ps13)$sample_type
  fung.rich.NoBee$temp_treat <- sample_data(ps13)$temp_treat
  fung.rich.NoBee$micro_treat <- sample_data(ps13)$micro_treat
  fung.rich.NoBee$combo_treat <- sample_data(ps13)$combo_treat
  fung.rich.NoBee$graft_stage <- sample_data(ps13)$graft_stage

# Plot richness and alpha diversity  
  phyloseq::plot_richness(ps13, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Examine the effect of temperature on Shannon diversity
  mod14 <- stats::aov(Shannon ~ temp_treat, data = fung.rich.NoBee)
  stats::anova(mod14)
  
# Post-hoc test
  TukeyHSD(mod14)
  
# Save p-values  
  stats.mod14 <- tibble::tribble(
                                 ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                   "WN",     "CN",    0.036,      "*"
                               )
  stats.mod14
  
# Examine the effect of temperature on Simpson diversity
  mod15 <- aov(Simpson ~ temp_treat, data = fung.rich.NoBee)
  stats::anova(mod15)
  
# Examine the effect of temperature on observed richness
  mod16 <- aov(Simpson ~ temp_treat, data = fung.rich.NoBee)
  stats::anova(mod16)
  
# Reorder x-axis
  fung.rich.NoBee$combo_treat <- factor(fung.rich.NoBee$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# Boxplot of Shannon index
  OsmiaCC.Shannon.fung.NoBee <- ggplot(fung.rich.NoBee, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                    values = climate.colors,
                                                    labels = climate.labs) +
                                  labs(title = fung.title) + 
                                  xlab("Treatment") +
                                  ylab("Shannon index") +
                                  ylim(0, 4) +
                                  ggpubr::stat_pvalue_manual(stats.mod14,
                                                             label = "p.adj.signif",
                                                             y.position = 3,
                                                             tip.length = 0.01)
  OsmiaCC.Shannon.fung.NoBee
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.fung.NoBee <- ggplot(fung.rich.NoBee, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = fung.title) + 
                                  xlab("Treatment") +
                                  ylab("Simpson index") +
                                  ylim(0, 1.0)
  OsmiaCC.Simpson.fung.NoBee
  
# Boxplot of Observed richness
  OsmiaCC.Observed.fung.NoBee <- ggplot(fung.rich.NoBee, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = fung.title) + 
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 25)
  OsmiaCC.Observed.fung.NoBee
  
# Provisions with bees  
  
# Estimate richness and alpha diversity
  fung.rich.bee <- phyloseq::estimate_richness(ps14, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  fung.rich.bee$sampleID <- sample_data(ps14)$sampleID
  fung.rich.bee$sample_type <- sample_data(ps14)$sample_type
  fung.rich.bee$temp_treat <- sample_data(ps14)$temp_treat
  fung.rich.bee$micro_treat <- sample_data(ps14)$micro_treat
  fung.rich.bee$combo_treat <- sample_data(ps14)$combo_treat
  fung.rich.bee$graft_stage <- sample_data(ps14)$graft_stage
  fung.rich.bee$sex <- sample_data(ps14)$sex
  
# Plot richness and alpha diversity
  phyloseq::plot_richness(ps14, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
    theme_bw() +
    xlab("")
  
# Examine the effects of sex, temperature, and microbiome treatment on Shannon diversity, with graft stage as a random effect
  mod17 <- nlme::lme(Shannon ~ sex + temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee)
  stats::anova(mod17)
  
# Post-hoc tests  
  emmeans(mod17, pairwise ~ sex, adjust = "tukey")
  emmeans(mod17, pairwise ~ micro_treat, adjust = "tukey")
  
# Examine the effects of sex, temperature, and microbiome treatment on Simpson diversity, with graft stage as a random effect
  mod18 <- nlme::lme(Simpson ~ sex + temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee)
  stats::anova(mod18)
  
# Post-hoc tests  
  emmeans(mod18, pairwise ~ sex, adjust = "tukey")
  emmeans(mod18, pairwise ~ micro_treat, adjust = "tukey")
  
# Examine the effects of sex, temperature, and microbiome treatment on observed richness, with graft stage as a random effect
  mod19 <- nlme::lme(Observed ~ sex + temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee)
  stats::anova(mod19)
  
# Post-hoc tests
  emmeans(mod19, pairwise ~ micro_treat, adjust = "tukey")
  
# Reorder x-axis
  fung.rich.bee$combo_treat <- factor(fung.rich.bee$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.fung.bee <- ggplot(fung.rich.bee, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                facet_grid(~ sex, 
                                           scale = "free", 
                                           space = "free") +
                                theme_bw() +
                                theme(plot.title = element_text(hjust = -0.12)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Treatment", 
                                                   values = climate.colors,
                                                   labels = climate.labs) +
                                labs(title = fung.title) +
                                xlab("Treatment") +
                                ylab("Shannon index") +
                                ylim(0, 3)
  OsmiaCC.Shannon.fung.bee
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.fung.bee <- ggplot(fung.rich.bee, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  facet_grid(~ sex, 
                                             scale = "free", 
                                             space = "free") +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = fung.title) + 
                                  xlab("Treatment") +
                                  ylab("Simpson index") +
                                  ylim(0, 1.0)
  OsmiaCC.Simpson.fung.bee
  
# Boxplot of Observed richness
  OsmiaCC.Observed.fung.bee <- ggplot(fung.rich.bee, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  facet_grid(~ sex, 
                                             scale = "free", 
                                             space = "free") +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment",
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = fung.title) +
                                  xlab("Treatment") +
                                  ylab("Observed richness") +
                                  ylim(0, 70)
  OsmiaCC.Observed.fung.bee
  
# Provisions with bees - males
  
# Subset data to include just males
  fung.rich.bee.M <- fung.rich.bee[fung.rich.bee$sex == "M", ]
  
# Examine the effect of temperature and microbiome treatment on Shannon diversity, with graft stage as a random effect
  mod20 <- nlme::lme(Shannon ~ temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee.M)
  stats::anova(mod20)
  
# Post-hoc test  
  emmeans(mod20, pairwise ~ micro_treat | temp_treat, adjust = "tukey")
  
# Save p-values  
  stats.mod20 <- tibble::tribble(
                                 ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                   "CS",     "CN",    0.0001,      "***",
                                   "AS",     "AN",    0.0001,      "***",
                                   "WS",     "WN",    0.0001,      "***"
                              )
  stats.mod20
  
# Examine the effect of temperature and microbiome treatment on Simpson diversity, with graft stage as a random effect
  mod21 <- nlme::lme(Simpson ~ temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee.M)
  stats::anova(mod21)
  
# Post-hoc test  
  emmeans(mod21, pairwise ~ micro_treat | temp_treat, adjust = "tukey")
  
# Save p-values  
  stats.mod21 <- tibble::tribble(
                                 ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                   "CS",     "CN",    0.0001,      "***",
                                   "AS",     "AN",    0.0001,      "***",
                                   "WS",     "WN",    0.0001,      "***"
                              )
  stats.mod21  
  
# Examine the effect of temperature and microbiome treatment on Observed richness, with graft stage as a random effect
  mod22 <- nlme::lme(Observed ~ temp_treat + micro_treat, random = ~1|graft_stage, data = fung.rich.bee.M)
  stats::anova(mod22)
  
# Post-hoc test  
  emmeans(mod22, pairwise ~ micro_treat | temp_treat, adjust = "tukey")
  
# Save p-values  
  stats.mod22 <- tibble::tribble(
                                 ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                   "CS",     "CN",    0.040,      "*",
                                   "AS",     "AN",    0.040,      "*",
                                   "WS",     "WN",    0.040,      "*"
                              )
  stats.mod22
  
# Reorder x-axis
  fung.rich.bee.M$combo_treat <- factor(fung.rich.bee.M$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.fung.bee.M <- ggplot(fung.rich.bee.M, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = fung.title) +
                                    xlab("Treatment") +
                                    ylab("Shannon index") +
                                    ylim(0, 4) + 
                                    ggpubr::stat_pvalue_manual(stats.mod20,
                                                               label = "p.adj.signif",
                                                               y.position = 3.5,
                                                               tip.length = 0.01)
  OsmiaCC.Shannon.fung.bee.M
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.fung.bee.M <- ggplot(fung.rich.bee.M, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = fung.title) + 
                                    xlab("Treatment") +
                                    ylab("Simpson index") +
                                    ylim(0, 1.0) +
                                    ggpubr::stat_pvalue_manual(stats.mod21,
                                                               label = "p.adj.signif",
                                                               y.position = 0.95,
                                                               tip.length = 0.01)
  OsmiaCC.Simpson.fung.bee.M
  
# Boxplot of Observed richness
  OsmiaCC.Observed.fung.bee.M <- ggplot(fung.rich.bee.M, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = fung.title) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 70) +
                                    ggpubr::stat_pvalue_manual(stats.mod22,
                                                               label = "p.adj.signif",
                                                               y.position = 65,
                                                               tip.length = 0.01)
  OsmiaCC.Observed.fung.bee.M
  
# Provisions with bees - females
  
# Subset data to include just females
  fung.rich.bee.F <- fung.rich.bee[fung.rich.bee$sex == "F", ]
  
# Examine the effect of temperature on Shannon diversity, with graft stage as a random effect
  mod23 <- nlme::lme(Shannon ~ temp_treat, random = ~1|graft_stage, data = fung.rich.bee.F)
  stats::anova(mod23)
  
# Examine the effect of temperature on Simpson diversity, with graft stage as a random effect
  mod24 <- nlme::lme(Simpson ~ temp_treat, random = ~1|graft_stage, data = fung.rich.bee.F)
  stats::anova(mod24)
  
# Examine the effect of temperature on Observed richness, with graft stage as a random effect
  mod25 <- nlme::lme(Observed ~ temp_treat, random = ~1|graft_stage, data = fung.rich.bee.F)
  stats::anova(mod25)
  
# Reorder x-axis
  fung.rich.bee.F$combo_treat <- factor(fung.rich.bee.F$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.fung.bee.F <- ggplot(fung.rich.bee.F, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = fung.title) +
                                  xlab("Treatment") +
                                  ylab("Shannon index") +
                                  ylim(0, 3.0)
  OsmiaCC.Shannon.fung.bee.F
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.fung.bee.F <- ggplot(fung.rich.bee.F, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = fung.title) + 
                                  xlab("Treatment") +
                                  ylab("Simpson index") +
                                  ylim(0, 1.0)
  OsmiaCC.Simpson.fung.bee.F
  
# Boxplot of Observed richness
  OsmiaCC.Observed.fung.bee.F <- ggplot(fung.rich.bee.F, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = fung.title) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 50)
  OsmiaCC.Observed.fung.bee.F
  
## Evenness ----
  
# Provisions with bees
  
# Extract ASV counts per sample
  otu <- phyloseq::otu_table(ps14)
  
# Calculate Shannon diversity index
  H <- vegan::diversity(otu, index = "shannon")
  
# Calculate observed richness
  S <- vegan::specnumber(otu)
  
# Calculate Pielou's evenness
  J <- H/log(S)
  
# Create df with diversity measures and metadata
  fung.evenness <- cbind(shannon = H, richness = S, pielou = J, sample_data(ps14))
  fung.evenness
  
# Remove samples with NaNs
  fung.evenness <- fung.evenness[complete.cases(fung.evenness), ]
  
# Examine the effects of sex, temperature, and microbiome treatment on evenness, with graft stage as a random effect
  mod26 <- nlme::lme(pielou ~ sex + temp_treat + micro_treat, random = ~1|graft_stage, data = fung.evenness)
  stats::anova(mod26)
  
# Plot
  OsmiaCC.Pielou.fung.bee <- ggplot(fung.evenness, aes(x = combo_treat, y = pielou, color = combo_treat)) +
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                facet_grid(~ sex, 
                                           scale = "free", 
                                           space = "free") +
                                theme_bw() +
                                theme(plot.title = element_text(hjust = -0.12)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                labs(title = fung.title) +
                                ylab("Pielou's Evenness") +
                                ylim(0, 1.0) +
                                xlab("") +
                                scale_color_manual(name = "Treatment",
                                                  values = climate.colors,
                                                  labels = climate.labs)
  OsmiaCC.Pielou.fung.bee

# Male bees
  
# Provisions with bees
  ps.fung.evenness.M <- subset_samples(ps14, sex == "M")
  
# Extract ASV counts per sample
  fung.otu.M <- phyloseq::otu_table(ps.fung.evenness.M)
  
# Calculate Shannon diversity index
  fung.H.M <- vegan::diversity(fung.otu.M, index = "shannon")
  
# Calculate observed richness
  fung.S.M <- vegan::specnumber(fung.otu.M)
  
# Calculate Pielou's evenness
  fung.J.M <- fung.H.M/log(fung.S.M)
  
# Create df with diversity measures and metadata
  fung.evenness.M <- cbind(shannon = fung.H.M, richness = fung.S.M, pielou = fung.J.M, sample_data(ps.fung.evenness.M))
  fung.evenness.M
  
# Remove samples with NaNs
  fung.evenness.M <- fung.evenness.M[complete.cases(fung.evenness.M), ]
  
# Examine the effects of temperature treatment and sex on evenness, with graft stage as a random effect
  mod27 <- nlme::lme(pielou ~ temp_treat + micro_treat, random = ~1|graft_stage, data = fung.evenness.M)
  stats::anova(mod27)
  
# Reorder x-axis
  fung.evenness.M$combo_treat <- factor(fung.evenness.M$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot
  OsmiaCC.Pielou.fung.bee.M <- ggplot(fung.evenness.M, aes(x = combo_treat, y = pielou, color = combo_treat)) +
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  labs(title = fung.title) +
                                  ylab("Pielou's Evenness") +
                                  ylim(0, 1.0) +
                                  xlab("Treatment") +
                                  scale_color_manual(name = "Treatment",
                                                     values = climate.colors,
                                                     labels = climate.labs)
  OsmiaCC.Pielou.fung.bee.M
  
## Beta diversity with relative abundance data ----
  
# Provisions with and without bees
  
# Calculate the relative abundance of each otu
  ps.prop.fung <- phyloseq::transform_sample_counts(ps15, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray <- phyloseq::distance(ps.prop.fung, method = "bray")
  
# Convert to df
  sample.fung <- data.frame(sample_data(ps15))
  
# Perform the PERMANOVA to test effects of temperature, microbiome, and sample type on fungal community composition
  fung.perm <- vegan::adonis2(fung.bray ~ temp_treat + micro_treat + sample_type, data = sample.fung)
  fung.perm
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.combo.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray, sample.fung$combo_treat, p.method = "BH")
  fung.perm.combo.BH
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.micro.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray, sample.fung$micro_treat, p.method = "BH")
  fung.perm.micro.BH
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.type.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray, sample.fung$sample_type, p.method = "BH")
  fung.perm.type.BH
  
# Provisions without bees
  
# Calculate the relative abundance of each otu  
  ps.prop.fung.NoBee <- phyloseq::transform_sample_counts(ps13, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.NoBee <- phyloseq::distance(ps.prop.fung.NoBee, method = "bray")
  
# Convert to df
  sample.fung.NoBee <- data.frame(sample_data(ps13))
  
# Perform the PERMANOVA to test effects of temperature on fungal community composition
  fung.perm.NoBee <- vegan::adonis2(fung.bray.NoBee ~ temp_treat, data = sample.fung.NoBee)
  fung.perm.NoBee
  
# Provisions with bees
  
# Calculate the relative abundance of each otu  
  ps.prop.fung.bee <- phyloseq::transform_sample_counts(ps14, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.bee <- phyloseq::distance(ps.prop.fung.bee, method = "bray")
  
# Convert to df
  sample.fung.bee <- data.frame(sample_data(ps14))
  
# Perform the PERMANOVA to test effects of temperature, microbiome, and sex on fungal community composition
  fung.perm.bee <- vegan::adonis2(fung.bray.bee ~ temp_treat + micro_treat + sex, data = sample.fung.bee)
  fung.perm.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.bee.micro.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray.bee, sample.fung.bee$micro_treat, p.method = "BH")
  fung.perm.bee.micro.BH

# Set permutations to deal with graft stage
  perm.relabund <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = sample.fung.bee$graft_stage,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of temperature, microbiome, and sex on fungal community composition, dealing with graft stage
  fung.perm.bee.graft <- vegan::adonis2(fung.bray.bee ~ temp_treat + micro_treat + sex, permutations = perm.relabund, data = sample.fung.bee)
  fung.perm.bee.graft
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.bee.micro.perm.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray.bee, sample.fung.bee$micro_treat, p.method = "BH")
  fung.perm.bee.micro.perm.BH
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# Provisions with and without bees
  
# Calculate the average distance of group members to the group centroid
  disp.fung.combo <- vegan::betadisper(fung.bray, sample.fung$combo_treat)
  disp.fung.combo
  
# Do any of the group dispersions differ?
  disp.fung.an.combo <- stats::anova(disp.fung.combo)
  disp.fung.an.combo
  
# Which group dispersions differ?
  disp.fung.tHSD.combo <- stats::TukeyHSD(disp.fung.combo)
  disp.fung.tHSD.combo
  
# Calculate the average distance of group members to the group centroid
  disp.fung.micro <- vegan::betadisper(fung.bray, sample.fung$micro_treat)
  disp.fung.micro
  
# Do any of the group dispersions differ?
  disp.fung.an.micro <- stats::anova(disp.fung.micro)
  disp.fung.an.micro

# Calculate the average distance of group members to the group centroid
  disp.fung.type <- vegan::betadisper(fung.bray, sample.fung$sample_type)
  disp.fung.type
  
# Do any of the group dispersions differ?
  disp.fung.an.type <- stats::anova(disp.fung.type)
  disp.fung.an.type

# Provisions without bees  
  
# Calculate the average distance of group members to the group centroid
  disp.fung.NoBee.temp <- vegan::betadisper(fung.bray.NoBee, sample.fung.NoBee$temp_treat)
  disp.fung.NoBee.temp
  
# Do any of the group dispersions differ?
  disp.fung.an.NoBee.temp <- stats::anova(disp.fung.NoBee.temp)
  disp.fung.an.NoBee.temp
  
# Provisions with bees
  
# Calculate the average distance of group members to the group centroid
  disp.fung.bee.combo <- vegan::betadisper(fung.bray.bee, sample.fung.bee$combo_treat)
  disp.fung.bee.combo
  
# Do any of the group dispersions differ?
  disp.fung.an.bee.combo <- stats::anova(disp.fung.bee.combo)
  disp.fung.an.bee.combo
  
# Which group dispersions differ?
  disp.fung.tHSD.combo <- stats::TukeyHSD(disp.fung.bee.combo)
  disp.fung.tHSD.combo
  
# Calculate the average distance of group members to the group centroid
  disp.fung.bee.temp <- vegan::betadisper(fung.bray.bee, sample.fung.bee$temp_treat)
  disp.fung.bee.temp
  
# Do any of the group dispersions differ?  
  disp.fung.an.bee.temp <- stats::anova(disp.fung.bee.temp)
  disp.fung.an.bee.temp
  
# Calculate the average distance of group members to the group centroid
  disp.fung.bee.micro <- vegan::betadisper(fung.bray.bee, sample.fung.bee$micro_treat)
  disp.fung.bee.micro
  
# Do any of the group dispersions differ?  
  disp.fung.an.bee.micro <- stats::anova(disp.fung.bee.micro)
  disp.fung.an.bee.micro
  
# Which group dispersions differ?
  disp.fung.tHSD.bee.micro <- stats::TukeyHSD(disp.fung.bee.micro)
  disp.fung.tHSD.bee.micro
  
# Calculate the average distance of group members to the group centroid
  disp.fung.bee.sex <- vegan::betadisper(fung.bray.bee, sample.fung.bee$sex)
  disp.fung.bee.sex
  
# Do any of the group dispersions differ?  
  disp.fung.an.bee.sex <- stats::anova(disp.fung.bee.sex)
  disp.fung.an.bee.sex
  
## Ordination with relative abundance data ----
  
# Provisions with and without bees
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung <- phyloseq::ordinate(ps.prop.fung, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung)$combo_treat <- factor(sample_data(ps.prop.fung)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung <- plot_ordination(ps.prop.fung, ord.pcoa.bray.fung, color = "combo_treat", shape = "sample_type") + 
                          theme_bw() +
                          theme(plot.title = element_text(hjust = -0.3)) +
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
                          theme(text = element_text(size = 16)) +
                          theme(legend.justification = "left", 
                                legend.title = element_text(size = 16, colour = "black"), 
                                legend.text = element_text(size = 14, colour = "black")) + 
                          geom_point(size = 3) +
                          scale_color_manual(values = climate.colors) +
                          labs(title = fung.title,
                               color = "Treatment",
                               shape = "Sample Type")
  OsmiaCC.PCoA.fung
  
# Provisions without bees
  
# PCoA using Bray-Curtis distance
  ord.pcoa.fung.NoBee <- phyloseq::ordinate(ps.prop.fung.NoBee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.NoBee)$combo_treat <- factor(sample_data(ps.prop.fung.NoBee)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.NoBee <- plot_ordination(ps.prop.fung.NoBee, ord.pcoa.fung.NoBee, color = "combo_treat") + 
                                theme_bw() +
                                theme(plot.title = element_text(hjust = -0.25)) +
                                theme(text = element_text(size = 16)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = fung.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.fung.NoBee
  
# Provisions with bees  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.fung.bee <- phyloseq::ordinate(ps.prop.fung.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.bee)$combo_treat <- factor(sample_data(ps.prop.fung.bee)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.bee <- plot_ordination(ps.prop.fung.bee, ord.pcoa.fung.bee, color = "combo_treat", shape = "sex") + 
                              theme_bw() +
                              theme(plot.title = element_text(hjust = -0.25)) +
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) +
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 16, colour = "black"), 
                                    legend.text = element_text(size = 14, colour = "black")) + 
                              geom_point(size = 3) +
                              scale_color_manual(values = climate.colors) +
                              labs(title = fung.title,
                                   color = "Treatment",
                                   shape = "Sex")
  OsmiaCC.PCoA.fung.bee

# Subset males
  ps.prop.fung.bee.M <- subset_samples(ps.prop.fung.bee, sex == "M")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung.bee.M <- phyloseq::ordinate(ps.prop.fung.bee.M, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.bee.M)$combo_treat <- factor(sample_data(ps.prop.fung.bee.M)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.bee.M <- plot_ordination(ps.prop.fung.bee.M, ord.pcoa.bray.fung.bee.M, color = "combo_treat") + 
                                theme_bw() +
                                theme(plot.title = element_text(hjust = -0.28)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(text = element_text(size = 16)) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = fung.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.fung.bee.M
  
# Subset females
  ps.prop.fung.bee.F <- subset_samples(ps.prop.fung.bee, sex == "F")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung.bee.F <- phyloseq::ordinate(ps.prop.fung.bee.F, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.bee.F)$combo_treat <- factor(sample_data(ps.prop.fung.bee.F)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.bee.F <- plot_ordination(ps.prop.fung.bee.F, ord.pcoa.bray.fung.bee.F, color = "combo_treat") + 
                                theme_bw() +
                                theme(plot.title = element_text(hjust = -0.26)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(text = element_text(size = 16)) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = fung.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.fung.bee.F
  
## Rarefaction ----
  
# Provisions with bees
  
# Produce rarefaction curves
  tab <- otu_table(ps14)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a tidy df
  rare.tidy.fungi <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  OsmiaCC.rare.fung.bee <- ggplot(rare.tidy.fungi, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(text = element_text(size = 10)) +
                            theme(plot.title = element_text(hjust = -0.11)) +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = fung.title) + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  OsmiaCC.rare.fung.bee

# Rarefy
  set.seed(1234)
  rareps.fung.bee <- rarefy_even_depth(ps14, sample.size = 25)
  
## Beta diversity with rarefied data ----  
  
# Provisions with bees
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.rare.bee <- phyloseq::distance(rareps.fung.bee, method = "bray")

# Convert to df
  sample.fung.rare.bee <- data.frame(sample_data(rareps.fung.bee))
  
# Perform the PERMANOVA to test effects of temperature, microbiome, and sex on fungal community composition 
  fung.perm.rare.bee <- vegan::adonis2(fung.bray.rare.bee ~ temp_treat + micro_treat + sex, data = sample.fung.rare.bee)
  fung.perm.rare.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.rare.bee.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray.rare.bee, sample.fung.rare.bee$micro_treat, p.method = "BH")
  fung.perm.rare.bee.BH
  
# Set permutations to deal with graft stage
  perm.rare.bee <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = sample.fung.rare.bee$graft_stage,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of temperature, microbiome, and sex on fungal community composition, using graft stage as a random effect
  fung.perm.rare.bee.graft <- vegan::adonis2(fung.bray.rare.bee ~ temp_treat + micro_treat + sex, permutations = perm.rare.bee, data = sample.fung.rare.bee)
  fung.perm.rare.bee.graft
  
# Follow up with pairwise comparisons - which sample types differ?
  fung.perm.rare.bee.graft.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray.rare.bee, sample.fung.rare.bee$micro_treat, p.method = "BH")
  fung.perm.rare.bee.graft.BH

## Test for homogeneity of multivariate dispersion with rarefied data ----
  
# Provisions with bees  
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare.bee.combo <- vegan::betadisper(fung.bray.rare.bee, sample.fung.rare.bee$combo_treat)
  disp.fung.rare.bee.combo
  
# Do any of the group dispersions differ?  
  disp.fung.an.rare.bee.combo <- stats::anova(disp.fung.rare.bee.combo)
  disp.fung.an.rare.bee.combo
  
# Which group dispersions differ?
  disp.fung.tHSD.rare.bee.combo <- stats::TukeyHSD(disp.fung.rare.bee.combo)
  disp.fung.tHSD.rare.bee.combo
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare.bee.temp <- vegan::betadisper(fung.bray.rare.bee, sample.fung.rare.bee$temp_treat)
  disp.fung.rare.bee.temp
  
# Do any of the group dispersions differ?  
  disp.fung.an.rare.bee.temp <- stats::anova(disp.fung.rare.bee.temp)
  disp.fung.an.rare.bee.temp
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare.bee.micro <- vegan::betadisper(fung.bray.rare.bee, sample.fung.rare.bee$micro_treat)
  disp.fung.rare.bee.micro
  
# Do any of the group dispersions differ?  
  disp.fung.an.rare.bee.micro <- stats::anova(disp.fung.rare.bee.micro)
  disp.fung.an.rare.bee.micro
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare.bee.sex <- vegan::betadisper(fung.bray.rare.bee, sample.fung.rare.bee$sex)
  disp.fung.rare.bee.sex
  
# Do any of the group dispersions differ?  
  disp.fung.an.rare.bee.sex <- stats::anova(disp.fung.rare.bee.sex)
  disp.fung.an.rare.bee.sex

## Ordination with rarefaction data ----

# Provisions with bees  
  
# Calculate the relative abundance of each otu  
  ps.prop.fung.rare.bee <- phyloseq::transform_sample_counts(rareps.fung.bee, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung.rare.bee <- phyloseq::ordinate(ps.prop.fung.rare.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.rare.bee)$combo_treat <- factor(sample_data(ps.prop.fung.rare.bee)$combo_treat, levels = c("CN", "AN", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.rare.bee <- plot_ordination(ps.prop.fung.rare.bee, ord.pcoa.bray.fung.rare.bee, color = "combo_treat", shape = "sex") + 
                                    theme_bw() +
                                    theme(plot.title = element_text(hjust = -0.25)) +
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 14, colour = "black")) + 
                                    geom_point(size = 3) +
                                    scale_color_manual(values = climate.colors) +
                                    labs(title = fung.title,
                                         color = "Treatment",
                                         shape = "Sex")
  OsmiaCC.PCoA.fung.rare.bee
  
# Subset males
  ps.prop.fung.rare.bee.M <- subset_samples(ps.prop.fung.rare.bee, sex == "M")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung.rare.bee.M <- phyloseq::ordinate(ps.prop.fung.rare.bee.M, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.rare.bee.M)$combo_treat <- factor(sample_data(ps.prop.fung.rare.bee.M)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.rare.bee.M <- plot_ordination(ps.prop.fung.rare.bee.M, ord.pcoa.bray.fung.rare.bee.M, color = "combo_treat") + 
                                      theme_bw() +
                                      theme(plot.title = element_text(hjust = -0.25)) +
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) + 
                                      geom_point(size = 3) +
                                      scale_color_manual(values = climate.colors) +
                                      labs(title = fung.title,
                                           color = "Treatment")
  OsmiaCC.PCoA.fung.rare.bee.M
  
# Subset females
  ps.prop.fung.rare.bee.F <- subset_samples(ps.prop.fung.rare.bee, sex == "F")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.fung.rare.bee.F <- phyloseq::ordinate(ps.prop.fung.rare.bee.F, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.rare.bee.F)$combo_treat <- factor(sample_data(ps.prop.fung.rare.bee.F)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.fung.rare.bee.F <- plot_ordination(ps.prop.fung.rare.bee.F, ord.pcoa.bray.fung.rare.bee.F, color = "combo_treat") + 
                                      theme_bw() +
                                      theme(plot.title = element_text(hjust = -0.26)) +
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) + 
                                      geom_point(size = 3) +
                                      scale_color_manual(values = climate.colors) +
                                      labs(title = fung.title,
                                           color = "Treatment")
  OsmiaCC.PCoA.fung.rare.bee.F
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  #Okabe.Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  #okabe.ext <- unikn::usecol(Okabe.Ito, n = 56)
  #colors <- sample(okabe.ext)
  
# Control provisions
  
# Agglomerate taxa by Genus
  y11 <- phyloseq::tax_glom(ps12, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y12 <- phyloseq::transform_sample_counts(y11, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y13 <- phyloseq::psmelt(y12)
  
# Ensure Genus is a chr
  y13$Genus <- as.character(y13$Genus)
  
# Group Genera with less that 1% abundance and rename
  y13$Genus[y13$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor
  y13$Genus <- as.factor(y13$Genus)
  
# Plot Genus by sample type 
  OsmiaCC.gen.fung.controls <- ggplot(data = y13, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                  geom_bar(stat = "identity", position = "fill") + 
                                  scale_fill_manual(values = colors) +
                                  theme(legend.position = "right",
                                        plot.title = element_text(hjust = -2.0)) +
                                  ylab("Relative abundance") + 
                                  ylim(0, 1.0) +
                                  xlab("Sample") +
                                  theme_bw() + 
                                  theme(text = element_text(size = 16)) +
                                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 18, colour = "black"), 
                                        legend.text = element_text(size = 16, colour = "black")) + 
                                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                  theme(panel.spacing.x = unit(0.1, "lines")) +
                                  guides(fill = guide_legend(ncol = 2)) +
                                  labs(fill = "Genera",
                                       title = fung.title)
  OsmiaCC.gen.fung.controls
  
# Provisions without bees
  
# Agglomerate taxa by Genus
  y14 <- phyloseq::tax_glom(ps13, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y15 <- phyloseq::transform_sample_counts(y14, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y16 <- phyloseq::psmelt(y15)
  
# Ensure Genus is a chr
  y16$Genus <- as.character(y16$Genus)
  
# Group Genera with less that 1% abundance and rename
  y16$Genus[y16$Abundance < 0.01] <- "Genus < 1% abund."
  
# Ensure Genus is a factor
  y16$Genus <- as.factor(y16$Genus)
  
# Reorder x-axis  
  y16$combo_treat <- factor(y16$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot treatment by Genus
  OsmiaCC.gen.fung.NoBee <- ggplot(data = y16, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                geom_bar(stat = "identity", position = "fill") + 
                                scale_fill_manual(values = colors) +
                                facet_grid(~ combo_treat, 
                                           scale = "free", 
                                           space = "free") +
                                theme(legend.position = "right",
                                      plot.title = element_text(hjust = -2.0)) +
                                ylab("Relative abundance") + 
                                ylim(0, 1.0) +
                                xlab("Treatment") +
                                theme_bw() + 
                                theme(text = element_text(size = 16)) +
                                theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank()) + 
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 18, colour = "black"), 
                                      legend.text = element_text(size = 16, colour = "black")) + 
                                guides(fill = guide_legend(ncol = 2)) +
                                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                theme(panel.spacing.x = unit(0.1, "lines")) +
                                labs(fill = "Genera",
                                     title = fung.title)
  OsmiaCC.gen.fung.NoBee
  
# Provisions with bees
  
# Agglomerate taxa by Genus
  y17 <- phyloseq::tax_glom(ps14, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y18 <- phyloseq::transform_sample_counts(y17, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y19 <- phyloseq::psmelt(y18)
  
# Ensure Genus is a chr
  y19$Genus <- as.character(y19$Genus)
  
# Group Genera with less that 1% abundance and rename
  y19$Genus[y19$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor
  y19$Genus <- as.factor(y19$Genus)
  
# Reorder x-axis  
  y19$combo_treat <- factor(y19$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot Genus by treatment
  OsmiaCC.gen.fung.bee <- ggplot(data = y19, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity", position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right",
                                    plot.title = element_text(hjust = -2.0)) +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Treatment") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 18, colour = "black"), 
                                    legend.text = element_text(size = 16, colour = "black")) + 
                              guides(fill = guide_legend(ncol = 3)) +
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              theme(panel.spacing.x = unit(0.1, "lines")) +
                              labs(fill = "Genera",
                                   title = fung.title)
  OsmiaCC.gen.fung.bee
  
# Subset data by sex
  y19.M <- y19[y19$sex == "M", ]
  y19.F <- y19[y19$sex == "F", ]
  
# Plot Genus for each sample - males
  OsmiaCC.gen.fung.bee.M <- ggplot(data = y19.M, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity",
                                       position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right",
                                    plot.title = element_text(hjust = -2.0)) +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Sample ID") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 18, colour = "black"), 
                                    legend.text = element_text(size = 16, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              theme(panel.spacing.x = unit(0.1, "lines")) +
                              guides(fill = guide_legend(ncol = 3)) +
                              labs(fill = "Genera",
                                   title = fung.title)
  OsmiaCC.gen.fung.bee.M
  
# Plot Genus for each sample - females  
  OsmiaCC.gen.fung.bee.F <- ggplot(data = y19.F, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity",
                                       position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right",
                                    plot.title = element_text(hjust = -2.0)) +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Sample ID") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 18, colour = "black"), 
                                    legend.text = element_text(size = 16, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              theme(panel.spacing.x = unit(0.1, "lines")) +
                              guides(fill = guide_legend(ncol = 3)) +
                              labs(fill = "Genera",
                                   title = fung.title)
  OsmiaCC.gen.fung.bee.F
  
# Top 15 Genera in male bees
  
# Provisions with bees
  y20 <- subset_samples(ps14, sex == "M")
  
# Agglomerate taxa by Genus
  y21 <- phyloseq::tax_glom(y20, taxrank = 'Genus')
  
# Identify the top 15 genera
  top15.fung.gen <- microbiome::top_taxa(y21, n = 15)
  
# Remove taxa that are not in the top 15
  ps.top15.fung.gen <- phyloseq::prune_taxa(top15.fung.gen, y21)
  
# Remove samples with 0 reads from the top 15 genera
  ps.top15.fung.gen <- phyloseq::prune_samples(sample_sums(ps.top15.fung.gen) != 0, ps.top15.fung.gen)
  
# Transform counts to relative abundances
  ps.top15.fung.gen.trans <- phyloseq::transform_sample_counts(ps.top15.fung.gen, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  ps.top15.fung.gen.trans.melt <- phyloseq::psmelt(ps.top15.fung.gen.trans)
  
# Order samples on x-axis
  ps.top15.fung.gen.trans.melt$combo_treat <- factor(ps.top15.fung.gen.trans.melt$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot top 15 genera for each sample
  OsmiaCC.15gen.relabund.fung <- ggplot(data = ps.top15.fung.gen.trans.melt, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ combo_treat, 
                                               scale = "free", 
                                               space = "free") +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") +
                                    ylim(0, 1.0) +
                                    scale_x_discrete(expand = c(0, 1)) +
                                    xlab("Sample") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) +
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 12, colour = "black"),
                                          strip.text = element_text(size = 14)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    theme(panel.spacing.x = unit(0.1, "lines")) +
                                    guides(fill = guide_legend(ncol = 1)) +
                                    labs(fill = "Genera",
                                         title = "B")
  OsmiaCC.15gen.relabund.fung

## Differential abundance with rarefied data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html

# Provisions with bees  
  
# Convert from a phyloseq to a deseq obj
  desq.obj.fung.rare.bee <- phyloseq::phyloseq_to_deseq2(ps14, ~ combo_treat)
  
# Calculate the geometric mean and remove rows with NA
  gm.mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means 
  geoMeans <- apply(counts(desq.obj.fung.rare.bee), 1, gm.mean)
  
# Estimate size factors
  desq.dds.fung.rare.bee <- estimateSizeFactors(desq.obj.fung.rare.bee, geoMeans = geoMeans)
  
# Fit a local regression
  desq.dds.fung.rare.bee <- DESeq(desq.dds.fung.rare.bee, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# CS vs CN
  
# Extract results from differential abundance table for CS vs CN
  CS.CN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "CS", "CN"))
  
# Order differential abundances by their padj value
  CS.CN.fung.rare.bee <- CS.CN.fung.rare.bee[order(CS.CN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  CS.CN.fung.rare.bee.p05 <- CS.CN.fung.rare.bee[(CS.CN.fung.rare.bee$padj < alpha & !is.na(CS.CN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  CS.CN.fung.rare.bee.p05
  
# AN vs CN
  
# Extract results from differential abundance table for AN vs CN
  AN.CN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "AN", "CN"))
  
# Order differential abundances by their padj value
  AN.CN.fung.rare.bee <- AN.CN.fung.rare.bee[order(AN.CN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AN.CN.fung.rare.bee.p05 <- AN.CN.fung.rare.bee[(AN.CN.fung.rare.bee$padj < alpha & !is.na(AN.CN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  AN.CN.fung.rare.bee.p05
  
# WN vs CN
  
# Extract results from differential abundance table for WN vs CN
  WN.CN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "WN", "CN"))
  
# Order differential abundances by their padj value
  WN.CN.fung.rare.bee <- WN.CN.fung.rare.bee[order(WN.CN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN.CN.fung.rare.bee.p05 <- WN.CN.fung.rare.bee[(WN.CN.fung.rare.bee$padj < alpha & !is.na(WN.CN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  WN.CN.fung.rare.bee.p05  
  
# AS vs AN
  
# Extract results from differential abundance table for AS vs AN
  AS.AN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "AS", "AN"))
  
# Order differential abundances by their padj value
  AS.AN.fung.rare.bee <- AS.AN.fung.rare.bee[order(AS.AN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AS.AN.fung.rare.bee.p05 <- AS.AN.fung.rare.bee[(AS.AN.fung.rare.bee$padj < alpha & !is.na(AS.AN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  AS.AN.fung.rare.bee.p05
  
# WN vs AN
  
# Extract results from differential abundance table for WN vs AN
  WN.AN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "WN", "AN"))
  
# Order differential abundances by their padj value
  WN.AN.fung.rare.bee <- WN.AN.fung.rare.bee[order(WN.AN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN.AN.fung.rare.bee.p05 <- WN.AN.fung.rare.bee[(WN.AN.fung.rare.bee$padj < alpha & !is.na(WN.AN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  WN.AN.fung.rare.bee.p05
  
# WS vs WN
  
# Extract results from differential abundance table for WS vs WN
  WS.WN.fung.rare.bee <- DESeq2::results(desq.dds.fung.rare.bee, contrast = c("combo_treat", "WS", "WN"))
  
# Order differential abundances by their padj value
  WS.WN.fung.rare.bee <- WS.WN.fung.rare.bee[order(WS.WN.fung.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WS.WN.fung.rare.bee.p05 <- WS.WN.fung.rare.bee[(WS.WN.fung.rare.bee$padj < alpha & !is.na(WS.WN.fung.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  WS.WN.fung.rare.bee.p05
  
