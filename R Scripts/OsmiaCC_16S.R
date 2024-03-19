##### Project: Osmia Climate Change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of 16S rRNA gene data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(stringr) # Version 1.5.1
  library(ggplot2) # Version 3.4.3
  library(plotrix) # Version 3.8-4
  library(phyloseq) # Version 1.44.0
  library(vegan) # Version 2.6-4
  library(RVAideMemoire) # Version 0.9-83-7
  library(magrittr) # Version 2.0.3
  library(tidyverse) # Version 1.2.0
  library(decontam) # Version 1.20.0
  library(microbiome) # Version 1.22.0
  library(nlme) # Version 3.1-164
  library(emmeans) # Version 1.10.0
  library(ggpubr) # Version 0.6.0
  library(car) # Version 3.1-2
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(DESeq2) # Version 1.40.2

# Italicize title
  bact.title <- expression(paste("(", italic("a"), ")"))

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

## Import data
  seqtab.nochim <- readRDS("OsmiaCC_seqs16Sall.rds")
  taxa <- readRDS("OsmiaCC_taxa16Sall.rds")
  meta16S.CC <- read.csv("OsmiaCC_master - 16S_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples.out <- stringr::str_sort(samples.out, numeric = TRUE)
  samples <- data.frame(meta16S.CC)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  temp_treat <- samples$temp_treat
  micro_treat <- samples$micro_treat
  combo_treat <- samples$combo_treat
  sample_or_control <- samples$sample_or_control
  sex <- samples$sex
  graft_stage <- samples$graft_stage
  DNA_conc <- samples$DNA_conc
  sample.info <- data.frame(extractionID = extractionID, 
                            sample_type = sample_type, 
                            sampleID = sampleID,  
                            temp_treat = temp_treat, 
                            micro_treat = micro_treat, 
                            combo_treat = combo_treat, 
                            sample_or_control = sample_or_control,
                            sex = sex,
                            graft_stage = graft_stage,
                            DNA_conc = DNA_conc)
  rownames(sample.info) <- samples.out
  
# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(sample.info), tax_table(taxa))
  ps1
  
# Display total number of reads, mean, and se in phyloseq obj before processing
  sum(sample_sums(ps1))
  mean(sample_sums(ps1))
  print(plotrix::std.error(sample_sums(ps1)))
  
## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
  
# Create df with LibrarySize for each sample
  df <- as.data.frame(sample_data(ps1))
  df$LibrarySize <- sample_sums(ps1)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  
# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_or_control)) + 
    geom_point()
  
# Determine which ASVs are contaminants based on frequency of DNA in negative controls
  contamdf.freq <- decontam::isContaminant(ps1, conc = DNA_conc, method = "frequency", threshold = 0.1)
  table(contamdf.freq$contaminant)
  head(which(contamdf.freq$contaminant))
  
# Determine which ASVs are contaminants based on frequency of DNA in negative controls with a higher threshold
  contamdf.freq05 <- decontam::isContaminant(ps1, conc = DNA_conc, method = "frequency", threshold = 0.5)
  table(contamdf.freq05$contaminant)
  head(which(contamdf.freq05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps1)$is.neg <- sample_data(ps1)$sample_or_control == "control"
  contamdf.prev <- decontam::isContaminant(ps1, method = "prevalence", neg = "is.neg", threshold = 0.1)
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls with a higher threshold
  contamdf.prev05 <- decontam::isContaminant(ps1, method = "prevalence", neg = "is.neg", threshold = 0.5)
  table(contamdf.prev05$contaminant)
  head(which(contamdf.prev05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls
  contamdf.comb <- decontam::isContaminant(ps1, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.1)
  table(contamdf.comb$contaminant)
  head(which(contamdf.comb$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls with a higher threshold  
  contamdf.comb05 <- decontam::isContaminant(ps1, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.5)
  table(contamdf.comb05$contaminant)
  head(which(contamdf.comb05$contaminant))
  
# Make phyloseq object of presence-absence in negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "control", ps1)
  
# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))
  
# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "sample", ps1)
  
# Calculate taxa abundance in samples from sample counts
  ps.pos.presence <- phyloseq::transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))
  
# Make data.frame of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contam.comb = contamdf.comb$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.comb)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")
  
# Make a new phyloseq object without contaminant taxa  
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.prev$contaminant, ps1)
  ps.noncontam
  
# Remove control samples used for identifying contaminants
  ps.sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps.sub
  
# Remove DNA from mitochondria & chloroplast
  ps2 <- ps.sub %>%
    phyloseq::subset_taxa(Kingdom == "Bacteria" &
                            Family  != "mitochondria" &
                              Class   != "Chloroplast"
    )
  
# Remove DNA from Eukarya, Eukaryota & Streptophyta
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Kingdom != "Eukarya" &
                          Kingdom != "Eukaryota" &
                            Family != "Streptophyta")
  
# Remove DNA from Archaea
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Kingdom != "Archaea")
  
# Remove samples without any reads
  ps3 <- phyloseq::prune_samples(sample_sums(ps2) != 0, ps2)
  ps3
  
# Transform counts to relative abundances
  ps3.relabund <- phyloseq::transform_sample_counts(ps3, function(x) x/sum(x))
  
# Save taxonomy, raw reads, and relative abundance data
  write.csv(tax_table(ps3), "OsmiaCC_16Staxa_all.csv")
  write.csv(otu_table(ps3), "OsmiaCC_16Sotu_all.csv")
  write.csv(otu_table(ps3.relabund), "OsmiaCC_16S_relabund.csv")
  
# Subset provisions collected before and after homogenization (controls)
  ps4 <- phyloseq::subset_samples(ps3, sample_type == "initial provision")
  ps4
  
# Provisions without bees  
  
# Subset data
  ps5 <- phyloseq::subset_samples(ps3, sample_type == "provision w/o bee")
  ps5
  
# Display total number of reads and means per sample in phyloseq obj after processing
  sum(sample_sums(ps5))
  mean(sample_sums(ps5))
  print(plotrix::std.error(sample_sums(ps5)))
  
# Calculate the reads per sample
  reads.sample.NoBee <- microbiome::readcount(ps5)
  head(reads.sample.NoBee)
  
# Add reads per sample to meta data
  sample_data(ps5)$reads.sample.NoBee <- reads.sample.NoBee
  
# Save sample metadata
  meta.NoBee <- sample_data(ps5)
  
# How many samples for each developmental stage?  
  meta.NoBee %>%
    group_by(sample_type, combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.NoBee),
              se = sd(reads.sample.NoBee)/sqrt(N),
              max = max(reads.sample.NoBee),
              min = min(reads.sample.NoBee))
  
# Add Seq to each taxa name
  taxa_names(ps5) <- paste0("Seq", seq(ntaxa(ps5)))
  
# Create a df containing the number of reads per OTU
  read.sums.df.NoBee <- data.frame(nreads = sort(taxa_sums(ps5), TRUE), 
                                                 sorted = 1:ntaxa(ps5), 
                                                 type = "OTUs")
  
# Add a column containing the number of reads per sample
  read.sums.df.NoBee <- rbind(read.sums.df.NoBee, data.frame(nreads = sort(sample_sums(ps5), TRUE), 
                                                             sorted = 1:nsamples(ps5), 
                                                             type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(read.sums.df.NoBee, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") + 
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")  
  
# Provisions with bees  
  
# Subset data
  ps6 <- phyloseq::subset_samples(ps3, sample_type == "final provision")
  ps6
  
# Display total number of reads and means per sample in phyloseq obj after processing
  sum(sample_sums(ps6))
  mean(sample_sums(ps6))
  print(plotrix::std.error(sample_sums(ps6)))
  
# Calculate the reads per sample
  reads.sample.bee <- microbiome::readcount(ps6)
  head(reads.sample.bee)
  
# Add reads per sample to meta data
  sample_data(ps6)$reads.sample.bee <- reads.sample.bee
  
# Save sample metadata
  meta.bee <- sample_data(ps6)

# Read stats for provisions with bees - males & females
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
  
# Read stats for provisions with bees - females
  meta.bee.F %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads.sample.bee),
              se = sd(reads.sample.bee)/sqrt(N),
              max = max(reads.sample.bee),
              min = min(reads.sample.bee))
  
# Add Seq to each taxa name
  taxa_names(ps6) <- paste0("Seq", seq(ntaxa(ps6)))
  
# Create a df containing the number of reads per OTU
  read.sums.df.bee <- data.frame(nreads = sort(taxa_sums(ps6), TRUE), 
                                 sorted = 1:ntaxa(ps6), 
                                 type = "OTUs")
  
# Add a column containing the number of reads per sample
  read.sums.df.bee <- rbind(read.sums.df.bee, data.frame(nreads = sort(sample_sums(ps6), TRUE), 
                                                         sorted = 1:nsamples(ps6), 
                                                         type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(read.sums.df.bee, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") + 
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
# Provisions with both male and female bees  
  
# Subset data
  ps7 <- phyloseq::subset_samples(ps3, sample_type != "initial provision")
  ps7
  
## Richness and alpha diversity ----  
  
# Provisions without bees  
  
# Estimate richness and alpha diversity
  bact.rich.NoBee <- phyloseq::estimate_richness(ps5, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  bact.rich.NoBee$sampleID <- sample_data(ps5)$sampleID
  bact.rich.NoBee$sample_type <- sample_data(ps5)$sample_type
  bact.rich.NoBee$temp_treat <- sample_data(ps5)$temp_treat
  bact.rich.NoBee$micro_treat <- sample_data(ps5)$micro_treat
  bact.rich.NoBee$combo_treat <- sample_data(ps5)$combo_treat
  
# Plot richness and alpha diversity 
  phyloseq::plot_richness(ps5, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Remove samples with 0 species richness 
  bact.rich.NoBee[bact.rich.NoBee == 0] <- NA
  bact.rich.NoBee <- bact.rich.NoBee[complete.cases(bact.rich.NoBee), ]

# Examine the effect of temperature on Shannon diversity
  mod1 <- stats::aov(Shannon ~ temp_treat, data = bact.rich.NoBee)
  stats::anova(mod1)

# Examine the effectsof temperature on Simpson diversity
  mod2 <- stats::aov(Simpson ~ temp_treat, data = bact.rich.NoBee)
  stats::anova(mod2)

# Examine the effect of temperature on observed richness
  mod3 <- stats::aov(Observed ~ temp_treat, data = bact.rich.NoBee)
  stats::anova(mod3)
  
# Reorder x-axis
  bact.rich.NoBee$combo_treat <- factor(bact.rich.NoBee$combo_treat, levels = c("CS", "CN","AS", "AN","WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.bact.NoBee <- ggplot(bact.rich.NoBee, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Shannon index") +
                                    ylim(0, 4)
  OsmiaCC.Shannon.bact.NoBee
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.bact.NoBee <- ggplot(bact.rich.NoBee, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) + 
                                    xlab("Treatment") +
                                    ylab("Simpson index") +
                                    ylim(0, 1.0)
  OsmiaCC.Simpson.bact.NoBee
  
# Boxplot of Observed richness
  OsmiaCC.Observed.bact.NoBee <- ggplot(bact.rich.NoBee, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 25)
  OsmiaCC.Observed.bact.NoBee

# Provisions with bees    
  
# Estimate Shannon, Simpson & observed richness
  bact.rich.bee <- phyloseq::estimate_richness(ps6, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bact.rich.bee$sampleID <- sample_data(ps6)$sampleID
  bact.rich.bee$sample_type <- sample_data(ps6)$sample_type
  bact.rich.bee$temp_treat <- sample_data(ps6)$temp_treat
  bact.rich.bee$micro_treat <- sample_data(ps6)$micro_treat
  bact.rich.bee$combo_treat <- sample_data(ps6)$combo_treat
  bact.rich.bee$graft_stage <- sample_data(ps6)$graft_stage
  bact.rich.bee$sex <- sample_data(ps6)$sex
  
# Plot richness and alpha diversity  
  phyloseq::plot_richness(ps6, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Remove samples with 0 species richness 
  bact.rich.bee[bact.rich.bee == 0] <- NA
  bact.rich.bee <- bact.rich.bee[complete.cases(bact.rich.bee), ]
  
# Examine the effects of temperature and sex on Shannon diversity, with graft stage as a random effect
  mod4 <- nlme::lme(Shannon ~ temp_treat + sex, random = ~1|graft_stage, data = bact.rich.bee)
  summary(mod4)
  stats::anova(mod4)
  
# Examine the effects of temperature and sex on Simpson diversity, with graft stage as a random effect
  mod5 <- nlme::lme(Simpson ~ temp_treat + sex, random = ~1|graft_stage, data = bact.rich.bee)
  summary(mod5)
  stats::anova(mod5)
  
# Examine the effects of temperature and sex on observed richness, with graft stage as a random effect
  mod6 <- nlme::lme(Observed ~ temp_treat + sex, random = ~1|graft_stage, data = bact.rich.bee)
  summary(mod6)
  stats::anova(mod6)
  
# Post-hoc test  
  emmeans(mod6, pairwise ~ sex + temp_treat, adjust = "tukey")
  
# Reorder x-axis
  bact.rich.bee$combo_treat <- factor(bact.rich.bee$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# Boxplot of Shannon index
  OsmiaCC.Shannon.bact.bee <- ggplot(bact.rich.bee, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(legend.position = "none",
                                        plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = bact.title) +
                                  xlab("Treatment") +
                                  ylab("Shannon index") +
                                  ylim(0, 4)
  OsmiaCC.Shannon.bact.bee
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.bact.bee <- ggplot(bact.rich.bee, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  facet_grid(~ sex, 
                                             scale = "free", 
                                             space = "free") +
                                  theme(legend.position = "none",
                                        plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate.colors,
                                                     labels = climate.labs) +
                                  labs(title = bact.title) + 
                                  xlab("Treatment") +
                                  ylab("Simpson index") +
                                  ylim(0, 1.0)
  OsmiaCC.Simpson.bact.bee
  
# Boxplot of Observed richness
  OsmiaCC.Observed.bact.bee <- ggplot(bact.rich.bee, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 70)
  OsmiaCC.Observed.bact.bee
  
# Provisions with bees - males
  
# Subset data to include just males
  bact.rich.bee.M <- bact.rich.bee[bact.rich.bee$sex == "M", ]
  
# Examine the effects of temperature on Shannon diversity, with graft stage as a random effect
  mod7 <- nlme::lme(Shannon ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.M)
  stats::anova(mod7)
  
# Examine the effects of temperature on Simpson diversity, with graft stage as a random effect
  mod8 <- nlme::lme(Simpson ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.M)
  stats::anova(mod8)
  
# Examine the effects of temperature on Observed richness, with graft stage as a random effect
  mod9 <- nlme::lme(Observed ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.M)
  stats::anova(mod9)
  
# Post-hoc test
  emmeans(mod9, pairwise ~ temp_treat, adjust = "tukey")
  
# Save p-values  
  stats.mod9 <- tibble::tribble(
                                ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                   "CN",     "WN",   0.0451,      "*"
                              )
  stats.mod9
  
# Reorder x-axis
  bact.rich.bee.M$combo_treat <- factor(bact.rich.bee.M$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.bact.bee.M <- ggplot(bact.rich.bee.M, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Shannon index") +
                                    ylim(0, 4)
  OsmiaCC.Shannon.bact.bee.M
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.bact.bee.M <- ggplot(bact.rich.bee.M, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) + 
                                    xlab("Treatment") +
                                    ylab("Simpson index") +
                                    ylim(0, 1.0)
  OsmiaCC.Simpson.bact.bee.M
  
# Boxplot of Observed richness
  OsmiaCC.Observed.bact.bee.M <- ggplot(bact.rich.bee.M, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 70) +
                                    labs(title = bact.title) +
                                    ggpubr::stat_pvalue_manual(stats.mod9,
                                                               label = "p.adj.signif",
                                                               y.position = 40,
                                                               tip.length = 0.01)
  OsmiaCC.Observed.bact.bee.M
  
# Provisions with bees - females
  
# Subset data to include just males
  bact.rich.bee.F <- bact.rich.bee[bact.rich.bee$sex == "F", ]
  
# Examine interactive effects of temperature on Shannon diversity, with graft stage as a random effect
  mod10 <- nlme::lme(Shannon ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.F)
  stats::anova(mod10)
  
# Examine interactive effects of temperature on Simpson diversity, with graft stage as a random effect
  mod11 <- nlme::lme(Simpson ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.F)
  stats::anova(mod11)
  
# Examine interactive effects of temperature on Observed richness, with graft stage as a random effect
  mod12 <- nlme::lme(Observed ~ temp_treat, random = ~1|graft_stage, data = bact.rich.bee.F)
  stats::anova(mod12)
  
# Reorder x-axis
  bact.rich.bee.F$combo_treat <- factor(bact.rich.bee.F$combo_treat, levels = c("CN","AN","WN"))
  
# Boxplot of Shannon index
  OsmiaCC.Shannon.bact.bee.F <- ggplot(bact.rich.bee.F, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Shannon index") +
                                    ylim(0, 3.0)
  OsmiaCC.Shannon.bact.bee.F
  
# Boxplot of Simpson index
  OsmiaCC.Simpson.bact.bee.F <- ggplot(bact.rich.bee.F, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment", 
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) + 
                                    xlab("Treatment") +
                                    ylab("Simpson index") +
                                    ylim(0, 1.0)
  OsmiaCC.Simpson.bact.bee.F
  
# Boxplot of Observed richness
  OsmiaCC.Observed.bact.bee.F <- ggplot(bact.rich.bee.F, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.12)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate.colors,
                                                       labels = climate.labs) +
                                    labs(title = bact.title) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ylim(0, 50)
  OsmiaCC.Observed.bact.bee.F

## Evenness ----
  
# Extract ASV counts per sample
  otu <- phyloseq::otu_table(ps6)
  
# Calculate Shannon diversity index
  H <- vegan::diversity(otu, index = "shannon")
  
# Calculate observed richness
  S <- vegan::specnumber(otu)
  
# Calculate Pielou's evenness
  J <- H/log(S)
  
# Create df with diversity measures and metadata
  bact.evenness <- cbind(shannon = H, richness = S, pielou = J, sample_data(ps6))
  bact.evenness
  
# Remove samples with NaNs
  bact.evenness <- bact.evenness[complete.cases(bact.evenness), ]

# Examine the effects of temperature treatment and sex on evenness, with graft stage as a random effect
  mod13 <- nlme::lme(pielou ~ temp_treat + sex, random = ~1|graft_stage, data = bact.evenness)
  summary(mod13)
  stats::anova(mod13)
  
# Plot
  OsmiaCC.Pielou.bact.bee <- ggplot(bact.evenness, aes(x = combo_treat, y = pielou, color = combo_treat)) +
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none",
                                      plot.title = element_text(hjust = -0.12)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                labs(title = bact.title) +
                                ylab("Pielou's Evenness") +
                                ylim(0, 1.0) +
                                xlab("Treatment") +
                                scale_color_manual(values = climate.colors,
                                                   labels = climate.labs)
  OsmiaCC.Pielou.bact.bee
  
# Male bees
  
# Provisions with bees
  ps.bact.evenness.M <- subset_samples(ps6, sex == "M")
  
# Extract ASV counts per sample
  bact.otu.M <- phyloseq::otu_table(ps.bact.evenness.M)
  
# Calculate Shannon diversity index
  bact.H.M <- vegan::diversity(bact.otu.M, index = "shannon")
  
# Calculate observed richness
  bact.S.M <- vegan::specnumber(bact.otu.M)
  
# Calculate Pielou's evenness
  bact.J.M <- bact.H.M/log(bact.S.M)
  
# Create df with diversity measures and metadata
  bact.evenness.M <- cbind(shannon = bact.H.M, richness = bact.S.M, pielou = bact.J.M, sample_data(ps.bact.evenness.M))
  bact.evenness.M
  
# Remove samples with NaNs
  bact.evenness.M <- bact.evenness.M[complete.cases(bact.evenness.M), ]
  
# Examine the effects of temperature treatment and sex on evenness, with graft stage as a random effect
  mod14 <- nlme::lme(pielou ~ temp_treat, random = ~1|graft_stage, data = bact.evenness.M)
  stats::anova(mod14)
  
# Reorder x-axis
  bact.evenness.M$combo_treat <- factor(bact.evenness.M$combo_treat, levels = c("CN", "AN","WN"))
  
# Plot
  OsmiaCC.Pielou.bact.bee.M <- ggplot(bact.evenness.M, aes(x = combo_treat, y = pielou, color = combo_treat)) +
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(legend.position = "none",
                                        plot.title = element_text(hjust = -0.12)) +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  labs(title = bact.title) +
                                  ylab("Pielou's Evenness") +
                                  ylim(0, 1.0) +
                                  xlab("Treatment") +
                                  scale_color_manual(values = climate.colors,
                                                     labels = climate.labs)
  OsmiaCC.Pielou.bact.bee.M
  
## Beta diversity with relative abundance data ----
  
# Provisions with and without bees
  
# Remove CS sample
  sample_data(ps7) <- sample_data(ps7)[!(sample_data(ps7)$combo_treat == "CS"), ]
  
# Calculate the relative abundance of each otu  
  ps.prop.bact <- phyloseq::transform_sample_counts(ps7, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray <- phyloseq::distance(ps.prop.bact, method = "bray")
  
# Convert to df
  sample.bact <- data.frame(sample_data(ps7))
  
# Perform the PERMANOVA to test effects of temperature and sample type on bacterial community composition
  bact.perm <- vegan::adonis2(bact.bray ~ temp_treat + sample_type, data = sample.bact)
  bact.perm
  
# Provisions without bees
  
# Calculate the relative abundance of each otu  
  ps.prop.bact.NoBee <- phyloseq::transform_sample_counts(ps5, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.NoBee <- phyloseq::distance(ps.prop.bact.NoBee, method = "bray")
  
# Convert to df
  sample.bact.NoBee <- data.frame(sample_data(ps5))
  
# Perform the PERMANOVA to test effects of treatments on bacterial community composition
  bact.perm.NoBee <- vegan::adonis2(bact.bray.NoBee ~ temp_treat, data = sample.bact.NoBee)
  bact.perm.NoBee
  
# Provisions with bees
  
# Remove CS sample
  sample_data(ps6) <- sample_data(ps6)[!(sample_data(ps6)$combo_treat == "CS"), ]

# Calculate the relative abundance of each otu  
  ps.prop.bact.bee <- phyloseq::transform_sample_counts(ps6, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.bee <- phyloseq::distance(ps.prop.bact.bee, method = "bray")
  
# Convert to df
  sample.bact.bee <- data.frame(sample_data(ps6))
  
# Perform the PERMANOVA to test effects of temperature and sex on bacterial community composition
  bact.perm.bee <- vegan::adonis2(bact.bray.bee ~ temp_treat + sex, data = sample.bact.bee)
  bact.perm.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.bee.combo.BH <- RVAideMemoire::pairwise.perm.manova(bact.bray.bee, sample.bact.bee$combo_treat, p.method = "BH")
  bact.perm.bee.combo.BH
  
# Set permutations to deal with graft stage
  perm.relabund <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = sample.bact.bee$graft_stage,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of temperature and sex on bacterial community composition, with graft stage as a random effect
  bact.perm.bee.graft <- vegan::adonis2(bact.bray.bee ~ temp_treat + sex, permutations = perm.relabund, data = sample.bact.bee)
  bact.perm.bee.graft
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# Provisions with and without bees

# Calculate the average distance of group members to the group centroid
  disp.bact.combo <- vegan::betadisper(bact.bray, sample.bact$combo_treat)
  disp.bact.combo
  
# Do any of the group dispersions differ?
  disp.bact.an.combo <- stats::anova(disp.bact.combo)
  disp.bact.an.combo
  
# Calculate the average distance of group members to the group centroid
  disp.bact.sex <- vegan::betadisper(bact.bray, sample.bact$sex)
  disp.bact.sex
  
# Do any of the group dispersions differ?
  disp.bact.an.sex <- stats::anova(disp.bact.sex)
  disp.bact.an.sex

# Provisions without bees  
  
# Calculate the average distance of group members to the group centroid
  disp.bact.NoBee.combo <- vegan::betadisper(bact.bray.NoBee, sample.bact.NoBee$combo_treat)
  disp.bact.NoBee.combo
  
# Do any of the group dispersions differ?
  disp.bact.an.NoBee.combo <- stats::anova(disp.bact.NoBee.combo)
  disp.bact.an.NoBee.combo
  
# Which group dispersions differ?
  disp.bact.tHSD.NoBee.combo <- stats::TukeyHSD(disp.bact.NoBee.combo)
  disp.bact.tHSD.NoBee.combo
  
# Provisions with bees  
  
# Calculate the average distance of group members to the group centroid
  disp.bact.bee.combo <- vegan::betadisper(bact.bray.bee, sample.bact.bee$combo_treat)
  disp.bact.bee.combo
  
# Do any of the group dispersions differ?
  disp.bact.an.bee.combo <- stats::anova(disp.bact.bee.combo)
  disp.bact.an.bee.combo

## Ordination with relative abundance data ----  
  
# Provisions with and without bees
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- phyloseq::ordinate(ps.prop.bact, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact)$combo_treat <- factor(sample_data(ps.prop.bact)$combo_treat, levels = c("CN", "AN", "WN"))

# Plot ordination
  OsmiaCC.PCoA.bact <- plot_ordination(ps.prop.bact, ord.pcoa.bray, color = "combo_treat", shape = "sample_type") + 
                          theme_bw() +
                          theme(legend.position = "none",
                                plot.title = element_text(hjust = -0.15)) +
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
                          theme(legend.justification = "left") + 
                          geom_point(size = 3) +
                          scale_color_manual(values = climate.colors) +
                          labs(title = bact.title,
                               color = "Treatment",
                               shape = "Sample Type")
  OsmiaCC.PCoA.bact
  
# Provisions without bees
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.NoBee <- phyloseq::ordinate(ps.prop.bact.NoBee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.NoBee)$combo_treat <- factor(sample_data(ps.prop.bact.NoBee)$combo_treat, levels = c("CN", "AN", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.NoBee <- plot_ordination(ps.prop.bact.NoBee, ord.pcoa.bray.NoBee, color = "combo_treat") + 
                                theme_bw() +
                                theme(legend.position = "none",
                                      plot.title = element_text(hjust = -0.25)) +
                                theme(text = element_text(size = 16)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = bact.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.bact.NoBee
  
# Provisions with bees  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bee <- phyloseq::ordinate(ps.prop.bact.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.bee)$combo_treat <- factor(sample_data(ps.prop.bact.bee)$combo_treat, levels = c("CN", "AN", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.bee <- plot_ordination(ps.prop.bact.bee, ord.pcoa.bray.bee, color = "combo_treat") + 
                                theme_bw() +
                                theme(legend.position = "none",
                                      plot.title = element_text(hjust = -0.25)) +
                                theme(text = element_text(size = 16)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = bact.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.bact.bee
  
# Subset males
  ps.prop.bact.bee.M <- subset_samples(ps.prop.bact.bee, sex == "M")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bact.bee.M <- phyloseq::ordinate(ps.prop.bact.bee.M, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.bee.M)$combo_treat <- factor(sample_data(ps.prop.bact.bee.M)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.bee.M <- plot_ordination(ps.prop.bact.bee.M, ord.pcoa.bray.bact.bee.M, color = "combo_treat") + 
                                theme_bw() +
                                theme(legend.position = "none",
                                      plot.title = element_text(hjust = -0.28)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(text = element_text(size = 16)) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = bact.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.bact.bee.M
  
# Subset females
  ps.prop.bact.bee.F <- subset_samples(ps.prop.bact.bee, sex == "F")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bact.bee.F <- phyloseq::ordinate(ps.prop.bact.bee.F, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.bee.F)$combo_treat <- factor(sample_data(ps.prop.bact.bee.F)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.bee.F <- plot_ordination(ps.prop.bact.bee.F, ord.pcoa.bray.bact.bee.F, color = "combo_treat") + 
                                theme_bw() +
                                theme(legend.position = "none",
                                      plot.title = element_text(hjust = -0.26)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(text = element_text(size = 16)) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate.colors) +
                                labs(title = bact.title,
                                     color = "Treatment")
  OsmiaCC.PCoA.bact.bee.F  
  
## Rarefaction ----
  
# Provisions with bees
  
# Produce rarefaction curves
  tab <- otu_table(ps6)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a tidy df
  rare.tidy.bact.bee <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  OsmiaCC.rare.bact.bee <- ggplot(rare.tidy.bact.bee, aes(x = Sample, y = Species, group = Site)) +
                              geom_line() +
                              theme_bw() +
                              theme(text = element_text(size = 10)) +
                              theme(plot.title = element_text(hjust = -0.11)) +
                              theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) +
                              labs(title = bact.title) + 
                              xlab("Number of reads") +
                              ylab("Number of species")
  OsmiaCC.rare.bact.bee
  
# Set seed and rarefy
  set.seed(1234)
  rareps.bact.bee <- phyloseq::rarefy_even_depth(ps6, sample.size = 15)
  
## Beta diversity with rarefied data ----  
  
# Provisions with bees  
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.rare.bee <- phyloseq::distance(rareps.bact.bee, method = "bray")
  
# Convert to df
  sample.bact.rare.bee <- data.frame(sample_data(rareps.bact.bee))
  
# Perform the PERMANOVA to test effects of temperature and sex on bacterial community composition  
  bact.perm.rare.bee <- vegan::adonis2(bact.bray.rare.bee ~ temp_treat + sex, data = sample.bact.rare.bee)
  bact.perm.rare.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.rare.bee.sex.BH <- RVAideMemoire::pairwise.perm.manova(bact.bray.rare.bee, sample.bact.rare.bee$sex, p.method = "BH")
  bact.perm.rare.bee.sex.BH
  
# Set permutations to deal with graft stage
  perm.rare.bee <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = sample.bact.rare.bee$graft_stage,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of temperature and sex on bacterial community composition, with graft stage as a random effect
  bact.perm.rare.bee.graft <- vegan::adonis2(bact.bray.rare.bee ~ temp_treat + sex, permutations = perm.rare.bee, data = sample.bact.rare.bee)
  bact.perm.rare.bee.graft

# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.rare.bee.sex.BH.graft <- RVAideMemoire::pairwise.perm.manova(bact.bray.rare.bee, sample.bact.rare.bee$sex, p.method = "BH")
  bact.perm.rare.bee.sex.BH.graft
  
## Test for homogeneity of multivariate dispersion with rarefied data ----
  
# Provisions with bees  
  
# Calculate the average distance of group members to the group centroid: combo_treat
  disp.bact.rare.bee.combo <- vegan::betadisper(bact.bray.rare.bee, sample.bact.rare.bee$combo_treat)
  disp.bact.rare.bee.combo
  
# Do any of the group dispersions differ?
  disp.bact.an.rare.bee.combo <- stats::anova(disp.bact.rare.bee.combo)
  disp.bact.an.rare.bee.combo
  
## Ordination with rarefied data ----
  
# Provisions with bees  
  
# Calculate the relative abundance of each otu  
  ps.prop.bact.rare.bee <- phyloseq::transform_sample_counts(rareps.bact.bee, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.rare.bee <- phyloseq::ordinate(ps.prop.bact.rare.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.rare.bee)$combo_treat <- factor(sample_data(ps.prop.bact.rare.bee)$combo_treat, levels = c("CN", "AN", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.rare.bee <- plot_ordination(ps.prop.bact.rare.bee, ord.pcoa.bray.rare.bee, color = "combo_treat", shape = "sex") + 
                                    theme_bw() +
                                    theme(legend.position = "none",
                                          plot.title = element_text(hjust = -0.25)) +
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 14, colour = "black")) + 
                                    geom_point(size = 3) +
                                    scale_color_manual(values = climate.colors) +
                                    labs(title = bact.title,
                                         color = "Treatment",
                                         shape = "Sex")
  OsmiaCC.PCoA.bact.rare.bee
  
# Subset males
  ps.prop.bact.rare.bee.M <- subset_samples(ps.prop.bact.rare.bee, sex == "M")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bact.rare.bee.M <- phyloseq::ordinate(ps.prop.bact.rare.bee.M, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.rare.bee.M)$combo_treat <- factor(sample_data(ps.prop.bact.rare.bee.M)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.rare.bee.M <- plot_ordination(ps.prop.bact.rare.bee.M, ord.pcoa.bray.bact.rare.bee.M, color = "combo_treat") + 
                                      theme_bw() +
                                      theme(legend.position = "none",
                                            plot.title = element_text(hjust = -0.28)) +
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) + 
                                      geom_point(size = 3) +
                                      scale_color_manual(values = climate.colors) +
                                      labs(title = bact.title,
                                           color = "Treatment")
  OsmiaCC.PCoA.bact.rare.bee.M
  
# Subset females
  ps.prop.bact.rare.bee.F <- subset_samples(ps.prop.bact.rare.bee, sex == "F")
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bact.rare.bee.F <- phyloseq::ordinate(ps.prop.bact.rare.bee.F, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.rare.bee.F)$combo_treat <- factor(sample_data(ps.prop.bact.rare.bee.F)$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot ordination
  OsmiaCC.PCoA.bact.rare.bee.F <- plot_ordination(ps.prop.bact.rare.bee.F, ord.pcoa.bray.bact.rare.bee.F, color = "combo_treat") + 
                                      theme_bw() +
                                      theme(legend.position = "none",
                                            plot.title = element_text(hjust = -0.26)) +
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) + 
                                      geom_point(size = 3) +
                                      scale_color_manual(values = climate.colors) +
                                      labs(title = bact.title,
                                           color = "Treatment")
  OsmiaCC.PCoA.bact.rare.bee.F
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  Okabe.Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe.ext <- unikn::usecol(Okabe.Ito, n = 63)
  colors <- sample(okabe.ext)
  
# Control provisions
  
# Agglomerate taxa by Genus
  y1 <- phyloseq::tax_glom(ps4, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y2 <- phyloseq::transform_sample_counts(y1, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y3 <- phyloseq::psmelt(y2)
  
# Ensure Genus is a chr
  y3$Genus <- as.character(y3$Genus)
  
# Group Genera with less that 1% abundance and rename
  y3$Genus[y3$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor
  y3$Genus <- as.factor(y3$Genus)
  
# Plot Genus by sample type 
  OsmiaCC.gen.bact.controls <- ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                  geom_bar(stat = "identity", position = "fill") + 
                                  scale_fill_manual(values = colors) +
                                  theme(legend.position = "right",
                                        plot.title = element_text(hjust = -2.0)) +
                                  ylab("Relative abundance") + 
                                  ylim(0, 1.0) +
                                  xlab("") +
                                  theme_bw() + 
                                  theme(text = element_text(size = 16)) +
                                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 18, colour = "black"), 
                                        legend.text = element_text(size = 16, colour = "black")) +
                                  theme(panel.spacing.x = unit(0.1, "lines")) +
                                  guides(fill = guide_legend(ncol = 1)) +
                                  labs(fill = "Genera",
                                       title = bact.title)
  OsmiaCC.gen.bact.controls

# Provisions without bees
  
# Agglomerate taxa by Genus
  y9 <- phyloseq::tax_glom(ps5, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y5 <- phyloseq::transform_sample_counts(y9, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y6 <- phyloseq::psmelt(y5)
  
# Ensure Genus is a chr
  y6$Genus <- as.character(y6$Genus)
  
# Group Genera with less that 1% abundance and rename
  y6$Genus[y6$Abundance < 0.01] <- "Genus < 1% abund."
  
# Ensure Genus is a factor
  y6$Genus <- as.factor(y6$Genus)
  
# Reorder x-axis  
  y6$combo_treat <- factor(y6$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot treatment by Genus
  OsmiaCC.gen.bact.NoBee <- ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity", position = "fill") + 
                              scale_fill_manual(values = colors) +
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              theme(legend.position = "right",
                                    plot.title = element_text(hjust = -2.0)) +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 18, colour = "black"), 
                                    legend.text = element_text(size = 16, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              theme(panel.spacing.x = unit(0.1, "lines")) +
                              guides(fill = guide_legend(ncol = 1)) +
                              labs(fill = "Genera",
                                   title = bact.title)
  OsmiaCC.gen.bact.NoBee
  
# Provisions with bees

# Agglomerate taxa by Genus
  y7 <- phyloseq::tax_glom(ps6, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y8 <- phyloseq::transform_sample_counts(y7, function(x) x/sum(x))
  
# Remove samples that have 0 reads
  y8 <- phyloseq::prune_samples(sample_sums(y8) != 0, y8)
  
# Convert to a ggplot2-friendly df
  y9 <- phyloseq::psmelt(y8)
  
# Ensure Genus is a chr
  y9$Genus <- as.character(y9$Genus)
  
# Group Genera with less that 1% abundance and rename
  y9$Genus[y9$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor
  y9$Genus <- as.factor(y9$Genus)
  
# Remove sample from CS
  y9 <- y9[y9$combo_treat != "CS", ]
  
# Reorder x-axis  
  y9$combo_treat <- factor(y9$combo_treat,levels = c("CN", "AN", "WN"))
  
# Plot Genus by treatment
  OsmiaCC.gen.bact.bee <- ggplot(data = y9, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity", position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right",
                                    plot.title = element_text(hjust = -2.0)) +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 18, colour = "black"), 
                                    legend.text = element_text(size = 16, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              theme(panel.spacing.x = unit(0.1, "lines")) +
                              guides(fill = guide_legend(ncol = 2)) +
                              labs(fill = "Genera",
                                   title = bact.title)
  OsmiaCC.gen.bact.bee
  
# Subset data by sex
  y9.M <- y9[y9$sex == "M", ]
  y9.F <- y9[y9$sex == "F", ]
  
# Plot Genus for each sample - males
  OsmiaCC.gen.bact.bee.M <- ggplot(data = y9.M, aes(x = sampleID, y = Abundance, fill = Genus)) + 
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
                              xlab("Sample") +
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
                                   title = bact.title)
  OsmiaCC.gen.bact.bee.M

# Plot Genus for each sample - females  
  OsmiaCC.gen.bact.bee.F <- ggplot(data = y9.F, aes(x = sampleID, y = Abundance, fill = Genus)) + 
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
                              xlab("Sample") +
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
                                   title = bact.title)
  OsmiaCC.gen.bact.bee.F
  
# Top 15 Genera in male bees
  
# Agglomerate taxa by Genus
  y10 <- phyloseq::tax_glom(ps6, taxrank = 'Genus')
  
# Identify the top 15 genera
  top15.bact.gen <- microbiome::top_taxa(y10, n = 15)
  
# Remove taxa that are not in the top 15
  ps.top15.bact.gen <- phyloseq::prune_taxa(top15.bact.gen, y10)
  
# Remove samples with 0 reads from the top 15 genera
  ps.top15.bact.gen <- phyloseq::prune_samples(sample_sums(ps.top15.bact.gen) != 0, ps.top15.bact.gen)
  
# Transform counts to relative abundances
  ps.top15.bact.gen.trans <- phyloseq::transform_sample_counts(ps.top15.bact.gen, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  ps.top15.bact.gen.trans.melt <- phyloseq::psmelt(ps.top15.bact.gen.trans)
  
# Order samples on x-axis
  ps.top15.bact.gen.trans.melt$combo_treat <- factor(ps.top15.bact.gen.trans.melt$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot top 15 genera for each sample
  OsmiaCC.15gen.relabund.bact <- ggplot(data = ps.top15.bact.gen.trans.melt, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ combo_treat, 
                                               scale = "free", 
                                               space = "free") +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    scale_x_discrete(expand = c(0, 1)) +
                                    xlab("") +
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
                                         title = bact.title)
  OsmiaCC.15gen.relabund.bact
  
## Differential abundance with rarefied data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# Provisions with bees  
  
# Convert from a phyloseq to a deseq obj
  desq.obj.bact.rare.bee <- phyloseq::phyloseq_to_deseq2(rareps.bact.bee, ~ combo_treat)
  
# Calculate the geometric mean and remove rows with NA
  gm.mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq.obj.bact.rare.bee), 1, gm.mean)
  
# Estimate size factors
  desq.dds.bact.rare.bee <- DESeq2::estimateSizeFactors(desq.obj.bact.rare.bee, geoMeans = geoMeans)
  
# Fit a local regression
  desq.dds.bact.rare.bee <- DESeq2::DESeq(desq.dds.bact.rare.bee, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05

# AN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  AN.CN.rare.bee <- DESeq2::results(desq.dds.bact.rare.bee, contrast = c("combo_treat", "AN", "CN"))
  
# Order differential abundances by their padj value
  AN.CN.rare.bee <- AN.CN.rare.bee[order(AN.CN.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AN.CN.rare.bee.p05 <- AN.CN.rare.bee[(AN.CN.rare.bee$padj < alpha & !is.na(AN.CN.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  AN.CN.rare.bee.p05
  
# WN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  WN.CN.rare.bee <- DESeq2::results(desq.dds.bact.rare.bee, contrast = c("combo_treat", "WN", "CN"))
  
# Order differential abundances by their padj value
  WN.CN.rare.bee <- WN.CN.rare.bee[order(WN.CN.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN.CN.rare.bee.p05 <- WN.CN.rare.bee[(WN.CN.rare.bee$padj < alpha & !is.na(WN.CN.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  WN.CN.rare.bee.p05
  
# WN vs AN
  
# Extract results from differential abundance table for initial vs final provision
  WN.AN.rare.bee <- DESeq2::results(desq.dds.bact.rare.bee, contrast = c("combo_treat", "WN", "AN"))
  
# Order differential abundances by their padj value
  WN.AN.rare.bee <- WN.AN.rare.bee[order(WN.AN.rare.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN.AN.rare.bee_p05 <- WN.AN.rare.bee[(WN.AN.rare.bee$padj < alpha & !is.na(WN.AN.rare.bee$padj)), ]
  
# Check to see if any padj is below alpha
  WN.AN.rare.bee_p05

# Arsenophonus ----
  
# Agglomerate taxa by Genus
  bact.gen <- phyloseq::tax_glom(ps6, taxrank = 'Genus')
  
# Transform to relative abundance  
  bact.gen <- phyloseq::transform_sample_counts(bact.gen, function(x) x/sum(x))
  
# Subset data to contain Arsenophonus only
  arsenophonus <- phyloseq::subset_taxa(bact.gen, Genus == "Arsenophonus")

# Save data
  write.csv(otu_table(arsenophonus), "OsmiaCC_Arsenophonus_otu.csv")

# Import data
  arseno <- read.csv("OsmiaCC_Arsenophonus.csv")
  
# Remove 0s
  arseno$rel_abund[arseno$rel_abund == 0] <- NA
  arseno <- arseno[complete.cases(arseno), ]
  
# Kruskal-Wallis test
  stats::kruskal.test(rel_abund ~ temp_treat, data = arseno)
  stats::kruskal.test(rel_abund ~ sex, data = arseno)
  
# Stats  
  arseno %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              mean = mean(rel_abund),
              se = sd(rel_abund)/sqrt(N))
  
# Reorder x-axis
  arseno$combo_treat <- factor(arseno$combo_treat, levels = c("CN", "AN","WN"))

# Plot
  arseno.rel.abund <- ggplot(arseno, aes(x = combo_treat, y = rel_abund, color = combo_treat)) + 
                          geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                          geom_jitter(size = 1, alpha = 0.9) +
                          theme_bw() +
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
                          scale_color_manual(name = "Treatment", 
                                             values = climate.colors) +
                          xlab("Treatment") +
                          ylab("Relative abundance") +
                          ylim(0, 1.0)
  arseno.rel.abund
  
