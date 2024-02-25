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
  library(car) # Version 3.1-2
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(DESeq2) # Version 1.40.2

# Import data
  seqtab.nochim <- readRDS("OsmiaCC_seqs16Sall.rds")
  taxa <- readRDS("OsmiaCC_taxa16Sall.rds")
  meta16S_CC <- read.csv("OsmiaCC_master - 16S_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples.out <- stringr::str_sort(samples.out, numeric = TRUE)
  samples <- data.frame(meta16S_CC)
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
  sampleinfo <- data.frame(extractionID = extractionID, 
                           sample_type = sample_type, 
                           sampleID = sampleID,  
                           temp_treat = temp_treat, 
                           micro_treat = micro_treat, 
                           combo_treat = combo_treat, 
                           sample_or_control = sample_or_control,
                           sex = sex,
                           graft_stage = graft_stage,
                           DNA_conc = DNA_conc)
  rownames(sampleinfo) <- samples.out
  
# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(sampleinfo), tax_table(taxa))
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
  ps_sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub
  
# Remove DNA from mitochondria & chloroplast
  ps2 <- ps_sub %>%
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
  
# Subset provisions collected before and after homogenization
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
  reads_sample <- microbiome::readcount(ps5)
  head(reads_sample)
  
# Add reads per sample to meta data
  sample_data(ps5)$reads_sample <- reads_sample
  
# Save sample metadata
  meta <- sample_data(ps5)
  
# How many samples for each developmental stage?  
  meta %>%
    group_by(sample_type, combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads_sample),
              se = sd(reads_sample)/sqrt(N),
              max = max(reads_sample),
              min = min(reads_sample))

# Save taxonomic and ASV counts
  write.csv(tax_table(ps5), "OsmiaCC_16Staxa_NoBees.csv")
  write.csv(otu_table(ps5), "OsmiaCC_16Sotu_NoBees.csv")
  
# Add Seq to each taxa name
  taxa_names(ps5) <- paste0("Seq", seq(ntaxa(ps5)))
  
# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps5), TRUE), 
                           sorted = 1:ntaxa(ps5), 
                           type = "OTUs")
  
# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps5), TRUE), 
                                             sorted = 1:nsamples(ps5), 
                                             type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
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
  reads_sample <- microbiome::readcount(ps6)
  head(reads_sample)
  
# Add reads per sample to meta data
  sample_data(ps6)$reads_sample <- reads_sample
  
# Save sample metadata
  meta <- sample_data(ps6)

# How many samples for each developmental stage?  
  meta %>%
    group_by(sample_type, combo_treat) %>%
    summarise(N = n(),
              mean = mean(reads_sample),
              se = sd(reads_sample)/sqrt(N),
              max = max(reads_sample),
              min = min(reads_sample))

# Save taxonomic and ASV counts
  write.csv(tax_table(ps6), "OsmiaCC_16Staxa_Bees.csv")
  write.csv(otu_table(ps6), "OsmiaCC_16Sotu_Bees.csv")
  
# Add Seq to each taxa name
  taxa_names(ps6) <- paste0("Seq", seq(ntaxa(ps6)))
  
# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps6), TRUE), 
                           sorted = 1:ntaxa(ps6), 
                           type = "OTUs")
  
# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps6), TRUE), 
                                             sorted = 1:nsamples(ps6), 
                                             type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") + 
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
## Richness and alpha diversity ----  
  
# Provisions without bees  
  
# Estimate Shannon, Simpson & observed richness
  bactrich <- phyloseq::estimate_richness(ps5, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bactrich$sampleID <- sample_data(ps5)$sampleID
  bactrich$sample_type <- sample_data(ps5)$sample_type
  bactrich$temp_treat <- sample_data(ps5)$temp_treat
  bactrich$micro_treat <- sample_data(ps5)$micro_treat
  bactrich$combo_treat <- sample_data(ps5)$combo_treat
  
# Plot Shannon index, Simpson index, & observed richness  
  phyloseq::plot_richness(ps5, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Remove samples with 0 species richness 
  bactrich[bactrich == 0] <- NA
  bactrich <- bactrich[complete.cases(bactrich), ]

# Examine interactive effects of temperature on Shannon diversity
  mod1 <- aov(Shannon ~ temp_treat, data = bactrich)
  stats::anova(mod1)

# Examine interactive effects of temperature on Simpson diversity
  mod2 <- aov(Simpson ~ temp_treat, data = bactrich)
  stats::anova(mod2)

# Examine interactive effects of temperature on observed richness
  mod3 <- aov(Observed ~ temp_treat, data = bactrich)
  stats::anova(mod3)
  
# Set color scheme
  climate_colors <- c("CS" = "#64B5F6",
                      "CN" = "#1565C0",
                      "AS" = "#9E9E9E",
                      "AN" = "#616161",
                      "WS" = "#E57373",
                      "WN" = "#C62828")
  
# Set labels
  climate_labs <- c("CS" = "Cool: Sterile",
                    "CN" = "Cool: Natural",
                    "AS" = "Ambient: Sterile",
                    "AN" = "Ambient: Natural",
                    "WS" = "Warm: Sterile",
                    "WN" = "Warm: Natural")
  
# Reorder x-axis
  bactrich$combo_treat <- factor(bactrich$combo_treat, levels = c("CS", "CN","AS", "AN","WS", "WN"))
  
# Boxplot of Shannon index
  OsmiaCC_Shannon_bact <- ggplot(bactrich, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              theme(legend.position = "none") +
                              theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) +
                              scale_color_manual(name = "Treatment", 
                                                 values = climate_colors,
                                                 labels = climate_labs) +
                              labs(title = "A") +
                              xlab("Treatment") +
                              ylab("Shannon index")
  OsmiaCC_Shannon_bact
  
# Boxplot of Simpson index
  OsmiaCC_Simpson_bact <- ggplot(bactrich, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              theme(legend.position = "none") +
                              theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) +
                              scale_color_manual(name = "Treatment", 
                                                 values = climate_colors,
                                                 labels = climate_labs) +
                              labs(title = "A") + 
                              xlab("Treatment") +
                              ylab("Simpson index")
  OsmiaCC_Simpson_bact
  
# Boxplot of Observed richness
  OsmiaCC_Observed_bact <- ggplot(bactrich, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              theme(legend.position = "none") +
                              theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) +
                              scale_color_manual(name = "Treatment",
                                                 values = climate_colors,
                                                 labels = climate_labs) +
                              xlab("Treatment") +
                              ylab("Observed richness") +
                              ggtitle("A")
  OsmiaCC_Observed_bact

# Provisions with bees    
  
# Estimate Shannon, Simpson & observed richness
  bactrich_bees <- phyloseq::estimate_richness(ps6, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bactrich_bees$sampleID <- sample_data(ps6)$sampleID
  bactrich_bees$sample_type <- sample_data(ps6)$sample_type
  bactrich_bees$temp_treat <- sample_data(ps6)$temp_treat
  bactrich_bees$micro_treat <- sample_data(ps6)$micro_treat
  bactrich_bees$combo_treat <- sample_data(ps6)$combo_treat
  bactrich_bees$graft_stage <- sample_data(ps6)$graft_stage
  bactrich_bees$sex <- sample_data(ps6)$sex
  
# Plot Shannon, Simpson & observed richness  
  phyloseq::plot_richness(ps6, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Remove samples with 0 species richness 
  bactrich_bees[bactrich_bees == 0] <- NA
  bactrich_bees <- bactrich_bees[complete.cases(bactrich_bees), ]
  
# Examine interactive effects of temperature on Shannon diversity
  mod4 <- nlme::lme(Shannon ~ temp_treat + sex, random = ~1|graft_stage, data = bactrich_bees)
  stats::anova(mod4)
  
# Examine interactive effects of temperature on Simpson diversity
  mod5 <- nlme::lme(Simpson ~ temp_treat + sex, random = ~1|graft_stage, data = bactrich_bees)
  stats::anova(mod5)
  
# Examine interactive effects of temperature on observed richness
  mod6 <- nlme::lme(Observed ~ temp_treat + sex, random = ~1|graft_stage, data = bactrich_bees)
  stats::anova(mod6)
  
  stats::shapiro.test(mod6$residuals)
  
  emmeans(mod6, pairwise ~ temp_treat, adjust = "tukey")
  emmeans(mod6, pairwise ~ sex, adjust = "tukey")
  
# Reorder x-axis
  bactrich_bees$combo_treat <- factor(bactrich_bees$combo_treat, levels = c("CN","AN","WN"))

# Boxplot of Shannon index
  OsmiaCC_Shannon_bact_bees <- ggplot(bactrich_bees, aes(x = combo_treat, y = Shannon, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(legend.position = "none") +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate_colors,
                                                     labels = climate_labs) +
                                  labs(title = "A") +
                                  xlab("Treatment") +
                                  ylab("Shannon index")
  OsmiaCC_Shannon_bact_bees
  
# Boxplot of Simpson index
  OsmiaCC_Simpson_bact_bees <- ggplot(bactrich_bees, aes(x = combo_treat, y = Simpson, color = combo_treat)) + 
                                  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                  geom_jitter(size = 1, alpha = 0.9) +
                                  theme_bw() +
                                  theme(legend.position = "none") +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  scale_color_manual(name = "Treatment", 
                                                     values = climate_colors,
                                                     labels = climate_labs) +
                                  labs(title = "A") + 
                                  xlab("Treatment") +
                                  ylab("Simpson index")
  OsmiaCC_Simpson_bact_bees
  
# Boxplot of Observed richness
  OsmiaCC_Observed_bact_bees <- ggplot(bactrich_bees, aes(x = combo_treat, y = Observed, color = combo_treat)) + 
                                    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                    geom_jitter(size = 1, alpha = 0.9) +
                                    theme_bw() +
                                    theme(legend.position = "none") +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                    scale_color_manual(name = "Treatment",
                                                       values = climate_colors,
                                                       labels = climate_labs) +
                                    xlab("Treatment") +
                                    ylab("Observed richness") +
                                    ggtitle("A")
  OsmiaCC_Observed_bact_bees
  
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
  bact_evenness <- cbind(shannon = H, richness = S, pielou = J, sample_data(ps6))
  bact_evenness
  
# Remove samples with NaNs
  bact_evenness <- bact_evenness[complete.cases(bact_evenness), ]
  
# Examine the effects of temperature treatment and sex on evenness using graft stage as a random effect
  mod7 <- nlme::lme(pielou ~ temp_treat + sex, random = ~1|graft_stage, data = bact_evenness)
  stats::anova(mod7)
  
# Plot
  OsmiaCC_Pielou_bact <- ggplot(bact_evenness, aes(x = combo_treat, y = pielou, color = combo_treat)) +
                            geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                            geom_jitter(size = 1, alpha = 0.9) +
                            theme_bw() +
                            theme(legend.position = "none") +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            ylab("Pielou's Evenness") +
                            xlab("Treatment") +
                            scale_color_manual(values = climate_colors,
                                               labels = climate_labs)
  OsmiaCC_Pielou_bact
  
## Beta diversity with relative abundance data ----
  
# Provisions without bees
  
# Calculate the relative abundance of each otu  
  ps.prop_bact_NoBees <- phyloseq::transform_sample_counts(ps5, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray_NoBees <- phyloseq::distance(ps.prop_bact_NoBees, method = "bray")
  
# Convert to data frame
  samplebact_NoBees <- data.frame(sample_data(ps5))
  
# Perform the PERMANOVA to test effects of treatments on bacterial community composition
  bact_perm_NoBees <- vegan::adonis2(bact_bray ~ temp_treat * micro_treat, data = samplebact)
  bact_perm_NoBees
  
# Provisions with bees
  
# Calculate the relative abundance of each otu  
  ps.prop_bact_bees <- phyloseq::transform_sample_counts(ps6, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray_bees <- phyloseq::distance(ps.prop_bact_bees, method = "bray")
  
# Convert to data frame
  samplebact_bees <- data.frame(sample_data(ps6))
  
# Perform the PERMANOVA to test effects of treatments on bacterial community composition
  bact_perm_bees <- vegan::adonis2(bact_bray_bees ~ temp_treat * micro_treat, data = samplebact_bees)
  bact_perm_bees
  
# Set permutations to deal with graft stage
  perm_relabund <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = samplebact_bees$graft_stage,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of treatments on bacterial community composition, dealing with graft stage
  bact_perm_graft <- vegan::adonis2(bact_bray_bees ~ temp_treat * micro_treat, permutations = perm_relabund, data = samplebact_bees)
  bact_perm_graft
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# Provisions without bees  
  
# Calculate the average distance of group members to the group centroid
  disp_bact_NoBees <- vegan::betadisper(bact_bray_NoBees, samplebact_NoBees$combo_treat)
  disp_bact_NoBees
  
# Do any of the group dispersions differ?
  disp_bact_an <- stats::anova(disp_bact_NoBees)
  disp_bact_an
  
# Which group dispersions differ?
  disp_bact_ttest_NoBees <- vegan::permutest(disp_bact_NoBees, 
                                             control = permControl(nperm = 999),
                                             pairwise = TRUE)
  disp_bact_ttest_NoBees
  
# Which group dispersions differ?
  disp_bact_tHSD_NoBees <- stats::TukeyHSD(disp_bact_NoBees)
  disp_bact_tHSD_NoBees
  
# Provisions with bees  
  
# Calculate the average distance of group members to the group centroid
  disp_bact <- vegan::betadisper(bact_bray, samplebact$combo_treat)
  disp_bact
  
# Do any of the group dispersions differ?
  disp_bact_an <- stats::anova(disp_bact)
  disp_bact_an
  
# Which group dispersions differ?
  disp_bact_ttest <- vegan::permutest(disp_bact, 
                                      control = permControl(nperm = 999),
                                      pairwise = TRUE)
  disp_bact_ttest
  
# Which group dispersions differ?
  disp_bact_tHSD <- stats::TukeyHSD(disp_bact)
  disp_bact_tHSD
  
# Calculate the average distance of group members to the group centroid: just temperature treatment
  disp_bact_temp <- vegan::betadisper(bact_bray, samplebact$temp_treat)
  disp_bact_temp
  
# Do any of the group dispersions differ?  
  disp_bact_temp_an <- stats::anova(disp_bact_temp)
  disp_bact_temp_an

## Ordination with relative abundance data ----  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_bees <- phyloseq::ordinate(ps.prop_bact_bees, method = "PCoA", distance = "bray")
  
# Plot ordination
  OsmiaCC_PCoA_bact_bees <- plot_ordination(ps.prop_bact_bees, ord.pcoa.bray_bees, color = "combo_treat", shape = "sample_type") + 
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(text = element_text(size = 16)) +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = climate_colors) +
                                labs(title = "A",
                                     color = "Treatment",
                                     shape = "Sample Type")
  OsmiaCC_PCoA_bact_bees
  
## Rarefaction ----
  
# Produce rarefaction curves
  tab <- otu_table(ps6)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a tidy df
  rare_tidy_bact <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  OsmiaCC_rare_bact <- ggplot(rare_tidy_bact, aes(x = Sample, y = Species, group = Site)) +
                          geom_line() +
                          theme_bw() +
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
                          labs(title = "A") + 
                          xlab("Number of reads") +
                          ylab("Number of species")
  OsmiaCC_rare_bact
  
# Set seed and rarefy
  #set.seed(1234)
  #rareps_bact <- phyloseq::rarefy_even_depth(ps6, sample.size = 10)
  
## Beta diversity with rarefied data ----  
  
# Create a distance matrix using Bray Curtis dissimilarity
  #bact_bray_rare <- phyloseq::distance(rareps_bact, method = "bray")
  
# Convert to data frame
  #samplebact_rare <- data.frame(sample_data(rareps_bact))
  
# Perform the PERMANOVA to test effects of treatment on bacterial community composition  
  #bact_perm_rare <- vegan::adonis2(bact_bray_rare ~ temp_treat, data = samplebact_rare)
  #bact_perm_rare
  
# Follow up with pairwise comparisons - which sample types differ?
  #bact_perm_rare_BH <- RVAideMemoire::pairwise.perm.manova(bact_bray_rare, samplebact_rare$combo_treat, p.method = "BH")
  #bact_perm_rare_BH
  
# Perform the PERMANOVA to test effects of treatment on bacterial community composition  
  #bact_perm_rare_graft <- vegan::adonis2(bact_bray_rare ~ graft_stage, data = samplebact_rare)
  #bact_perm_rare_graft
  
# Set permutations to deal with graft stage
  #perm_rare <- how(within = Within(type = "free"),
                   #plots = Plots(type = "none"),
                   #blocks = samplebact_rare$graft_stage,
                   #observed = FALSE,
                   #complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with graft stage
  #bact_perm_rare_graft <- vegan::adonis2(bact_bray_rare ~ temp_treat, permutations = perm_rare, data = samplebact_rare)
  #bact_perm_rare_graft
  
# Follow up with pairwise comparisons - which sample types differ?
  #bact_perm_rare_BH <- RVAideMemoire::pairwise.perm.manova(bact_bray_rare, samplebact_rare$combo_treat, p.method = "BH")
  #bact_perm_rare_BH
  
## Test for homogeneity of multivariate dispersion with rarefied data ----
  
# Calculate the average distance of group members to the group centroid: combo_treat
  #disp_bact_rare_combo <- vegan::betadisper(bact_bray_rare, samplebact_rare$combo_treat)
  #disp_bact_rare_combo
  
# Do any of the group dispersions differ?
  #disp_bact_an_rare <- stats::anova(disp_bact_rare_combo)
  #disp_bact_an_rare
  
# Which group dispersions differ?
  #disp_bact_ttest_rare_combo <- vegan::permutest(disp_bact_rare_combo, 
                                                 #control = permControl(nperm = 999),
                                                 #pairwise = TRUE)
  #disp_bact_ttest_rare_combo
  
# Which group dispersions differ?
  #disp_bact_tHSD_rare_combo <- stats::TukeyHSD(disp_bact_rare_combo)
  #disp_bact_tHSD_rare_combo
  
## Ordination with rarefied data ----
  
# Calculate the relative abundance of each otu  
  #ps.prop_rare <- phyloseq::transform_sample_counts(rareps_bact, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  #ord.pcoa.bray_rare <- phyloseq::ordinate(ps.prop_rare, method = "PCoA", distance = "bray")
  
# Plot ordination
  #OsmiaCC_PCoA_bact_rare <- plot_ordination(ps.prop_rare, ord.pcoa.bray_rare, color = "combo_treat", shape = "sample_type") + 
                                #theme_bw() +
                                #theme(legend.position = "none") +
                                #theme(text = element_text(size = 16)) +
                                #theme(legend.justification = "left", 
                                      #legend.title = element_text(size = 16, colour = "black"), 
                                      #legend.text = element_text(size = 14, colour = "black")) + 
                                #geom_point(size = 3) +
                                #scale_color_manual(values = c("#616161", "#9E9E9E", "#1565C0", "#64B5F6", "#C62828", "#E57373")) +
                                #labs(title = "A", color = "Treatment", shape = "Sample Type")
  #OsmiaCC_PCoA_bact_rare
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe_ext <- unikn::usecol(Okabe_Ito, n = 57)
  colors <- sample(okabe_ext)
  
# Provisions without bees
  
# Sort data by Family
  y15 <- phyloseq::tax_glom(ps5, taxrank = 'Family') # agglomerate taxa
  y16 <- phyloseq::transform_sample_counts(y15, function(x) x/sum(x))
  y16 <- phyloseq::psmelt(y16)
  y16$Family <- as.character(y16$Family)
  y16$Family[y16$Abundance < 0.01] <- "Family < 1% abund."
  y16$Family <- as.factor(y16$Family)
  head(y16)
  
# Plot treatment by Family
  ggplot(data = y16, aes(x = combo_treat, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 1)) +
    ggtitle("A")
  
# Sort data by Genus
  y17 <- phyloseq::tax_glom(ps5, taxrank = 'Genus') # agglomerate taxa
  y18 <- phyloseq::transform_sample_counts(y17, function(x) x/sum(x))
  y18 <- phyloseq::psmelt(y18)
  y18$Genus <- as.character(y18$Genus)
  y18$Genus[y18$Abundance < 0.01] <- "Genus < 1% abund."
  y18$Genus <- as.factor(y18$Genus)
  head(y18)
  
# Plot treatment by Genus
  ggplot(data = y18, aes(x = combo_treat, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 1)) +
    ggtitle("A")
  
# Provisions with bees  
  
# Sort data by Family
  y1 <- phyloseq::tax_glom(ps6, taxrank = 'Family') # agglomerate taxa
  y2 <- phyloseq::transform_sample_counts(y1, function(x) x/sum(x))
  y3 <- phyloseq::psmelt(y2)
  y3$Family <- as.character(y3$Family)
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  y3$Family <- as.factor(y3$Family)
  head(y3)
  
# Remove provision from CS
  y3 <- y3[y3$combo_treat != "CS", ]
  
# Save relative abundance data
  write.csv(y3, "OsmiaCC_Fam_bact_relabund.csv")
  
# Reorder x-axis  
  y3$combo_treat <- factor(y3$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot treatment by Family
  OsmiaCC_fam_relabund_bact <- ggplot(data = y3, aes(x = combo_treat, y = Abundance, fill = Family)) + 
                                  geom_bar(stat = "identity", position = "fill") + 
                                  scale_fill_manual(values = colors) + 
                                  theme(legend.position = "right") +
                                  ylab("Relative abundance") + 
                                  ylim(0, 1.0) +
                                  xlab("Treatment") +
                                  theme_bw() + 
                                  theme(text = element_text(size = 16)) +
                                  theme(panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank()) + 
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 16, colour = "black"), 
                                        legend.text = element_text(size = 14, colour = "black")) + 
                                  guides(fill = guide_legend(ncol = 2)) +
                                  ggtitle("A")
  OsmiaCC_fam_relabund_bact
  
# Plot Family for each sample
  ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) +
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample ID") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Bacteria")

# Sort data by Genus
  y4 <- phyloseq::tax_glom(ps6, taxrank = 'Genus') # agglomerate taxa
  y5 <- phyloseq::transform_sample_counts(y4, function(x) x/sum(x))
  y5 <- phyloseq::prune_samples(sample_sums(y5) != 0, y5)
  y6 <- phyloseq::psmelt(y5)
  y6$Genus <- as.character(y6$Genus)
  y6$Genus[y6$Abundance < 0.01] <- "Genera < 1% abund."
  y6$Genus <- as.factor(y6$Genus)
  head(y6)
  
# Save relative abundance data
  write.csv(y6, "OsmiaCC_Gen_bact_relabund.csv")
  
# Remove sample from CS  
  y6 <- y6[y6$combo_treat != "CS", ]
  
# Reorder x-axis  
  y6$combo_treat <- factor(y6$combo_treat,levels = c("CS", "CN", "AN", "WN"))
  
# Plot Genus by treatment
  OsmiaCC_gen_type_bact <- ggplot(data = y6, aes(x = combo_treat, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity", position = "fill") + 
                              facet_grid(~ sex, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right") +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Treatment") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 16, colour = "black"), 
                                    legend.text = element_text(size = 14, colour = "black")) + 
                              guides(fill = guide_legend(ncol = 2)) +
                              ggtitle("A")
  OsmiaCC_gen_type_bact
  
# Subset by sex
  y6_males <- y6[y6$sex == "M", ]
  y6_females <- y6[y6$sex == "F", ]
  
# Plot Genus for each sample - males
  OsmiaCC_gen_ID_bact_M <- ggplot(data = y6_males, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity",
                                       position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right") +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Sample ID") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 16, colour = "black"), 
                                    legend.text = element_text(size = 14, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              guides(fill = guide_legend(ncol = 2)) +
                              ggtitle("A")
  OsmiaCC_gen_ID_bact_M
  
  ggsave("OsmiaCC_16Sgenera_males.png", plot = OsmiaCC_gen_ID_bact_M, width = 15, height = 10, unit = "in")
  
# Plot Genus for each sample - females  
  OsmiaCC_gen_ID_bact_F <- ggplot(data = y6_females, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                              geom_bar(stat = "identity",
                                       position = "fill") + 
                              facet_grid(~ combo_treat, 
                                         scale = "free", 
                                         space = "free") +
                              scale_fill_manual(values = colors) +
                              theme(legend.position = "right") +
                              ylab("Relative abundance") + 
                              ylim(0, 1.0) +
                              xlab("Sample ID") +
                              theme_bw() + 
                              theme(text = element_text(size = 16)) +
                              theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) + 
                              theme(legend.justification = "left", 
                                    legend.title = element_text(size = 16, colour = "black"), 
                                    legend.text = element_text(size = 14, colour = "black")) + 
                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                              guides(fill = guide_legend(ncol = 2)) +
                              ggtitle("A")
  OsmiaCC_gen_ID_bact_F
  
  ggsave("OsmiaCC_16Sgenera_females.png", plot = OsmiaCC_gen_ID_bact_F, width = 15, height = 10, unit = "in")
  
# Abundances
  
# Agglomerate taxa
  
  bact_gen <- phyloseq::tax_glom(ps6, taxrank = 'Genus')
  bact_gen <- phyloseq::transform_sample_counts(bact_gen, function(x) x/sum(x))
  bact_gen <- phyloseq::prune_samples(sample_sums(bact_gen) != 0, bact_gen)
  write.csv(tax_table(bact_gen), "OsmiaCC_bact_gen_taxa.csv")
  write.csv(otu_table(bact_gen), "OsmiaCC_bact_gen_otu.csv")
  
  #bact_abund <- read.csv("OsmiaCC_bact_abundance.csv")
  
  
  
  
  
  
  
  
  
  
  
## Differential abundance with raw data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# Convert from a phyloseq to a deseq obj
  desq_obj_rare <- phyloseq::phyloseq_to_deseq2(rareps_bact, ~ combo_treat)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq_obj_rare), 1, gm_mean)
  
# Estimate size factors
  desq_dds_rare <- DESeq2::estimateSizeFactors(desq_obj_rare, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds_rare <- DESeq2::DESeq(desq_dds_rare, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# WN vs AN
  
# Extract results from differential abundance table for initial vs final provision
  WN_AN_rare <- DESeq2::results(desq_dds_rare, contrast = c("combo_treat", "WN", "AN"))
  
# Order differential abundances by their padj value
  WN_AN_rare <- WN_AN_rare[order(WN_AN_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_AN_rare_p05 <- WN_AN_rare[(WN_AN_rare$padj < alpha & !is.na(WN_AN_rare$padj)), ]
  
# Check to see if any padj is below alpha
  WN_AN_rare_p05
  
# AN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  AN_CN_rare <- DESeq2::results(desq_dds_rare, contrast = c("combo_treat", "AN", "CN"))
  
# Order differential abundances by their padj value
  AN_CN_rare <- AN_CN_rare[order(AN_CN_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AN_CN_rare_p05 <- AN_CN_rare[(AN_CN_rare$padj < alpha & !is.na(AN_CN_rare$padj)), ]
  
# Check to see if any padj is below alpha
  AN_CN_rare_p05
  
# WN vs CN

# Extract results from differential abundance table for initial vs final provision
  WN_CN_rare <- DESeq2::results(desq_dds_rare, contrast = c("combo_treat", "WN", "CN"))

# Order differential abundances by their padj value
  WN_CN_rare <- WN_CN_rare[order(WN_CN_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_CN_rare_p05 <- WN_CN_rare[(WN_CN_rare$padj < alpha & !is.na(WN_CN_rare$padj)), ]
  
# Check to see if any padj is below alpha
  WN_CN_rare_p05
  
