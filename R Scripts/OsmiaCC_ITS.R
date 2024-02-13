##### Project: Osmia Climate Change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of ITS data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(phyloseq) # Version 1.44.0
  library(vegan) # Version 2.6-4
  library(magrittr) # Version 2.0.3
  library(tidyverse) # Version 1.2.0
  library(decontam) # Version 1.20.0
  library(nlme) # Version 3.1-164
  library(emmeans) # Version 1.10.0
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(DESeq2) # Version 1.40.2

# Import data
  seqtab.nochim1 <- readRDS("OsmiaCC_seqsITS_run1.rds")
  taxa1 <- readRDS("OsmiaCC_taxaITS_run1.rds")
  metaITS_CC_run1 <- read.csv("OsmiaCC_master - ITS_run1.csv")
  
  seqtab.nochim2 <- readRDS("OsmiaCC_seqsITS_run2.rds")
  taxa2 <- readRDS("OsmiaCC_taxaITS_run2.rds")
  metaITS_CC_run2 <- read.csv("OsmiaCC_master - ITS_run2.csv")

## Create phyloseq objects for each ITS run ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim1)
  samples <- data.frame(metaITS_CC_run1)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  temp_treat <- samples$temp_treat
  micro_treat <- samples$micro_treat
  combo_treat <- samples$combo_treat
  sample_or_control <- samples$sample_or_control
  sex <- samples$sex
  graft_stage <- samples$graft_stage
  sampleinfo <- data.frame(extractionID = extractionID, 
                           sample_type = sample_type, 
                           sampleID = sampleID,  
                           temp_treat = temp_treat, 
                           micro_treat = micro_treat, 
                           combo_treat = combo_treat, 
                           sample_or_control = sample_or_control,
                           sex = sex,
                           graft_stage = graft_stage)
  rownames(sampleinfo) <- samples.out

# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows = FALSE), sample_data(sampleinfo), tax_table(taxa1))
  ps1
  
# Re-create your df
  samples.out <- rownames(seqtab.nochim2)
  samples <- data.frame(metaITS_CC_run2)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  temp_treat <- samples$temp_treat
  micro_treat <- samples$micro_treat
  combo_treat <- samples$combo_treat
  sample_or_control <- samples$sample_or_control
  sex <- samples$sex
  graft_stage <- samples$graft_stage
  sampleinfo <- data.frame(extractionID = extractionID,
                           sample_type = sample_type,
                           sampleID = sampleID,
                           temp_treat = temp_treat,
                           micro_treat = micro_treat,
                           combo_treat = combo_treat,
                           sample_or_control = sample_or_control,
                           sex = sex,
                           graft_stage = graft_stage)
  rownames(sampleinfo) <- samples.out
  
# Format your data to work with phyloseq
  ps2 <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows = FALSE), sample_data(sampleinfo), tax_table(taxa2))
  ps2
  
# Merge phyloseq objects
  ps3 <- merge_phyloseq(ps1, ps2)
  ps3
  
# Display total number of reads and means per sample in phyloseq obj before processing
  sum(sample_sums(ps3))
  mean(sample_sums(ps3))
  
## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Create df with LibrarySize for each sample 
  df <- as.data.frame(sample_data(ps3))
  df$LibrarySize <- sample_sums(ps3)
  df <- df[order(df$LibrarySize), ]
  df$Index <- seq(nrow(df))
  
# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_or_control)) + 
    geom_point()
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps3)$is.neg <- sample_data(ps3)$sample_or_control == "control"
  contamdf.prev <- decontam::isContaminant(ps3, method = "prevalence", neg = "is.neg", threshold = 0.1)
  
# How many contaminants are there?
  table(contamdf.prev$contaminant)
  
# Which ASVs are contaminants?
  head(which(contamdf.prev$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) higher than 0.5 in negative controls
  contamdf.prev05 <- decontam::isContaminant(ps3, method = "prevalence", neg = "is.neg", threshold = 0.5)
  table(contamdf.prev05$contaminant)
  
# Which ASVs are contaminants?
  head(which(contamdf.prev05$contaminant))
  
# Make phyloseq object of presence-absence in negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps3)$sample_or_control == "control", ps3)
  
# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))
  
# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps3)$sample_or_control == "sample", ps3)
  
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
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.prev05$contaminant, ps3)
  ps.noncontam
  
# Remove control samples used for identifying contaminants
  ps_sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub
  
# Remove samples without any reads
  ps4 <- phyloseq::prune_samples(sample_sums(ps_sub) != 0, ps_sub)
  ps4
  
# Remove samples from females
  ps4 <- phyloseq::prune_samples(sample_data(ps4)$sex != "F", ps4)
  ps4
  
# Display total number of reads and means per sample in phyloseq obj after processing
  sum(sample_sums(ps4))
  mean(sample_sums(ps4))
  
# Save sample metadata
  meta <- sample_data(ps4)
  
# How many samples for each developmental stage?  
  meta %>%
    group_by(sample_type, combo_treat) %>%
    summarise(N = n())
  
# Save taxonomic and ASV counts
  write.csv(tax_table(ps4), "OsmiaCC_ITStaxa.csv")
  write.csv(otu_table(ps4), "OsmiaCC_ITSotu.csv")
  
# Add Seq to each taxa name
  taxa_names(ps4) <- paste0("Seq", seq(ntaxa(ps4)))
  
# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps4), TRUE), 
                           sorted = 1:ntaxa(ps4),
                           type = "OTUs")
  
# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps4), TRUE), 
                                             sorted = 1:nsamples(ps4),
                                             type = "Samples"))
  
# Plot number of reads per ASV and sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") +
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
## Alpha diversity ----
  
# Calculate species richness
  fungrich <- phyloseq::estimate_richness(ps4, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  fungrich$sampleID <- sample_data(ps4)$sampleID
  fungrich$sample_type <- sample_data(ps4)$sample_type
  fungrich$temp_treat <- sample_data(ps4)$temp_treat
  fungrich$micro_treat <- sample_data(ps4)$micro_treat
  fungrich$combo_treat <- sample_data(ps4)$combo_treat
  fungrich$graft_stage <- sample_data(ps4)$graft_stage

# Plot species richness  
  phyloseq::plot_richness(ps4, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
                theme_bw() +
                xlab("")
  
# Remove samples with 0 species richness
  fungrich[fungrich == 0] <- NA
  fungrich <- fungrich[complete.cases(fungrich), ]
  
# Examine interactive effects of temperature and microbiome treatments on Shannon diversity
  mod4 <- nlme::lme(Shannon ~ temp_treat * micro_treat, random = ~1|graft_stage, data = fungrich)
  anova(mod4)
  
# Examine interactive effects of temperature and microbiome treatments on Simpson diversity
  mod5 <- nlme::lme(Simpson ~ temp_treat * micro_treat, random = ~1|graft_stage, data = fungrich)
  anova(mod5)
  
# Examine interactive effects of temperature and microbiome treatments on observed richness
  mod6 <- nlme::lme(Simpson ~ temp_treat * micro_treat, random = ~1|graft_stage, data = fungrich)
  anova(mod6)
  
# Reorder x-axis
  fungrich$combo_treat <- factor(fungrich$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# New names for facet_grid
  type_names <- c('final provision' = "provisions with bee",
                  'provision w/o bee' = "provisions without bee")

# Boxplot of Shannon index
  OsmiaCC_Shannon_fungi <- ggplot(fungrich, aes(x = combo_treat, y = Shannon, fill = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              facet_grid(~ sample_type,
                                         scale = "free",
                                         space = "free",
                                         labeller = as_labeller(type_names)) +
                              scale_fill_manual(name = "Treatment", 
                                                values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                                labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                              labs(title = "A") + 
                              xlab("Treatment") +
                              ylab("Shannon index")
  OsmiaCC_Shannon_fungi
  
# Boxplot of Simpson index
  OsmiaCC_Simpson_fungi <- ggplot(fungrich, aes(x = combo_treat, y = Simpson, fill = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              facet_grid(~ sample_type,
                                         scale = "free",
                                         space = "free",
                                         labeller = as_labeller(type_names)) +
                              scale_fill_manual(name = "Treatment", 
                                                values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                                labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                              labs(title = "B") + 
                              xlab("Treatment") +
                              ylab("Simpson index")
  OsmiaCC_Simpson_fungi
  
# Boxplot of Observed richness
  OsmiaCC_Observed_fungi <- ggplot(fungrich, aes(x = combo_treat, y = Observed, fill = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              facet_grid(~ sample_type,
                                         scale = "free",
                                         space = "free",
                                         labeller = as_labeller(type_names)) +
                              scale_fill_manual(name = "Treatment", 
                                                values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                                labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                              labs(title = "B") + 
                              xlab("Treatment") +
                              ylab("Observed richness")
  OsmiaCC_Observed_fungi
  
## Beta diversity with relative abundance data ----
  
# Calculate the relative abundance of each otu  
  ps.prop_fung <- phyloseq::transform_sample_counts(ps4, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung_bray <- phyloseq::distance(ps.prop_fung, method = "bray")
  
# Convert to data frame
  samplefung <- data.frame(sample_data(ps4))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  fung_perm <- vegan::adonis2(fung_bray ~ temp_treat * micro_treat, data = samplefung)
  fung_perm
  
# Follow up with pairwise comparisons - which sample types differ? by temperature treatment
  fung_perm_temp_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray, samplefung$temp_treat, p.method = "BH")
  fung_perm_temp_BH
  
# Follow up with pairwise comparisons - which sample types differ? by microbiome treatment  
  fung_perm_micro_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray, samplefung$micro_treat, p.method = "BH")
  fung_perm_micro_BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_relabund <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = samplefung$graft_stage,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  fung_perm <- vegan::adonis2(fung_bray ~ temp_treat * micro_treat, permutations = perm_relabund, data = samplefung)
  fung_perm
  
# Follow up with pairwise comparisons - which sample types differ? by temperature treatment
  fung_perm_temp_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray, samplefung$temp_treat, p.method = "BH")
  fung_perm_temp_BH
  
# Follow up with pairwise comparisons - which sample types differ? by microbiome treatment  
  fung_perm_micro_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray, samplefung$micro_treat, p.method = "BH")
  fung_perm_micro_BH
  
## Test for homogeneity of multivariate dispersion without rarefaction ----
  
# Calculate the average distance of group members to the group centroid
  disp_fung <- vegan::betadisper(fung_bray, samplefung$combo_treat)
  disp_fung
  
# Do any of the group dispersions differ?
  disp_fung_an <- anova(disp_fung)
  disp_fung_an
  
# Calculate the average distance of group members to the group centroid: just temperature treatment
  disp_fung_temp <- vegan::betadisper(fung_bray, samplefung$temp_treat)
  disp_fung_temp
  
# Do any of the group dispersions differ?  
  disp_fung_temp_an <- anova(disp_fung_temp)
  disp_fung_temp_an
  
# Calculate the average distance of group members to the group centroid: just microbiome treatment
  disp_fung_micro <- vegan::betadisper(fung_bray, samplefung$micro_treat)
  disp_fung_micro
  
# Do any of the group dispersions differ?  
  disp_fung_micro_an <- anova(disp_fung_micro)
  disp_fung_micro_an
  
# Which group dispersions differ?
  disp_fung_ttest <- vegan::permutest(disp_fung, 
                                      control = permControl(nperm = 999),
                                      pairwise = TRUE)
  disp_fung_ttest
  
# Which group dispersions differ?
  disp_fung_tHSD <- stats::TukeyHSD(disp_fung)
  disp_fung_tHSD  
  
## Ordination with relative abundance data ----
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_fung <- phyloseq::ordinate(ps.prop_fung, method = "PCoA", distance = "bray")
  
# Plot ordination
  OsmiaCC_PCoA_fung <- plot_ordination(ps.prop_fung, ord.pcoa.bray_fung, color = "combo_treat", shape = "sample_type") + 
                          theme_bw() +
                          theme(legend.position = "none") +
                          theme(text = element_text(size = 16)) +
                          theme(legend.justification = "left", 
                                legend.title = element_text(size = 16, colour = "black"), 
                                legend.text = element_text(size = 14, colour = "black")) + 
                          geom_point(size = 3) +
                          scale_color_manual(values = c("#616161", "#9E9E9E", "#1565C0", "#64B5F6", "#C62828", "#E57373")) +
                          labs(title = "B", color = "Treatment", shape = "Sample Type")
  OsmiaCC_PCoA_fung
  
## Rarefaction ----
  
# Produce rarefaction curves
  tab <- otu_table(ps4)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a tidy df
  rare_tidy_fungi <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  OsmiaCC_rare_fungi <- ggplot(rare_tidy_fungi, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "B") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  OsmiaCC_rare_fungi

# Rarefy
  set.seed(1234)
  fung_rareps <- rarefy_even_depth(ps4, sample.size = 11)
  
# Perform PERMANOVA to test effects of treatments on bacterial community composition
  fung_bray_rare <- phyloseq::distance(fung_rareps, method = "bray")
  
# Convert to data frame
  samplefung_rare <- data.frame(sample_data(fung_rareps))
  
# Perform the PERMANOVA to test effects of treatment on bacterial community composition 
  fung_perm_rare <- vegan::adonis2(fung_bray_rare ~ temp_treat * micro_treat, data = samplefung_rare)
  fung_perm_rare
  
# Follow up with pairwise comparisons - which sample types differ? microbiome treatment only
  fung_perm_rare_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray_rare, samplefung_rare$micro_treat, p.method = "BH")
  fung_perm_rare_BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_rare <- how(within = Within(type = "free"),
                   plots = Plots(type = "none"),
                   blocks = samplefung_rare$graft_stage,
                   observed = FALSE,
                   complete = FALSE)
  
# Perform the PERMANOVA to test effects of treatment on bacterial community composition 
  fung_perm_rare <- vegan::adonis2(fung_bray_rare ~ temp_treat * micro_treat, permutations = perm_rare, data = samplefung_rare)
  fung_perm_rare
  
# Follow up with pairwise comparisons - which sample types differ? microbiome treatment only
  fung_perm_rare_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray_rare, samplefung_rare$micro_treat, p.method = "BH")
  fung_perm_rare_BH

## Test for homogeneity of multivariate dispersion ----
  
# Calculate the average distance of group members to the group centroid
  disp_fung_rare <- vegan::betadisper(fung_bray_rare, samplefung_rare$combo_treat)
  disp_fung_rare
  
# Do any of the group dispersions differ?
  disp_fung_an_rare <- anova(disp_fung_rare)
  disp_fung_an_rare
  
# Calculate the average distance of group members to the group centroid: just temperature treatment
  disp_fung_temp_rare <- vegan::betadisper(fung_bray_rare, samplefung_rare$temp_treat)
  disp_fung_temp_rare
  
# Do any of the group dispersions differ?  
  disp_fung_temp_an_rare <- anova(disp_fung_temp_rare)
  disp_fung_temp_an_rare
  
# Calculate the average distance of group members to the group centroid: just microbiome treatment
  disp_fung_micro_rare <- vegan::betadisper(fung_bray_rare, samplefung_rare$micro_treat)
  disp_fung_micro_rare
  
# Do any of the group dispersions differ?  
  disp_fung_micro_an_rare <- anova(disp_fung_micro_rare)
  disp_fung_micro_an_rare
  
# Which group dispersions differ?
  disp_fung_ttest_rare <- vegan::permutest(disp_fung_rare, 
                                           control = permControl(nperm = 999),
                                           pairwise = TRUE)
  disp_fung_ttest_rare
  
# Which group dispersions differ?
  disp_fung_tHSD_rare <- stats::TukeyHSD(disp_fung_rare)
  disp_fung_tHSD_rare
  
## Ordination ----

# Calculate the relative abundance of each otu  
  ps.prop_rare <- phyloseq::transform_sample_counts(fung_rareps, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_rare <- phyloseq::ordinate(ps.prop_rare, method = "PCoA", distance = "bray")
  
# Plot ordination
  OsmiaCC_PCoA_fungi_rare <- plot_ordination(ps.prop_rare, ord.pcoa.bray_rare, color = "combo_treat", shape = "sample_type") + 
                                theme_bw() +
                                theme(text = element_text(size = 16)) +
                                theme(legend.justification = "left", 
                                      legend.title = element_text(size = 16, colour = "black"), 
                                      legend.text = element_text(size = 14, colour = "black")) + 
                                geom_point(size = 3) +
                                scale_color_manual(values = c( "#616161", "#9E9E9E", "#1565C0", "#64B5F6", "#C62828", "#E57373")) +
                                labs(title = "B", color = "Treatment", shape = "Sample Type")
  OsmiaCC_PCoA_fungi_rare
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe_ext <- unikn::usecol(Okabe_Ito, n = 41)
  colors <- sample(okabe_ext)
  
# Remove patterns in tax_table   
  tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))] <- gsub(tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))], pattern = "[a-z]__", replacement = "")
  
# Sort data by Family
  y7 <- phyloseq::tax_glom(fung_rareps, taxrank = 'Family') # agglomerate taxa
  y8 <- phyloseq::transform_sample_counts(y7, function(x) x/sum(x))
  y9 <- phyloseq::psmelt(y8)
  y9$Family <- as.character(y9$Family)
  y9$Family[y9$Abundance < 0.01] <- "Family < 1% abund."
  y9$Family <- as.factor(y9$Family)
  head(y9)
  
# Save relative abundance data
  write.csv(y9, "OsmiaCC_Fam_bact_relabund.csv")  
  
# Reorder x-axis  
  y9$combo_treat <- factor(y9$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# New names for facet_grid
  type_names <- c('final provision' = "provisions with bee",
                  'provision w/o bee' = "provisions without bee")
  
# Plot Family by treatment  
  OsmiaCC_fam_relabund_fungi <- ggplot(data = y9, aes(x = combo_treat, y = Abundance, fill = Family)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) + 
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free",
                                               labeller = as_labeller(type_names)) +
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
                                    ggtitle("Fungi")
  OsmiaCC_fam_relabund_fungi
  
# Plot Family for each sample
  ggplot(data = y9, aes(x = sampleID, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    facet_grid(~ combo_treat, 
               scale = "free", 
               space = "free") +
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Fungi")
  
# Sort data by Genus
  y10 <- phyloseq::tax_glom(fung_rareps, taxrank = 'Genus') # agglomerate taxa
  y11 <- phyloseq::transform_sample_counts(y10, function(x) x/sum(x))
  y12 <- phyloseq::psmelt(y11)
  y12$Genus <- as.character(y12$Genus)
  y12$Genus[y12$Abundance < 0.01] <- "Genera < 1% abund."
  y12$Genus <- as.factor(y12$Genus)
  head(y12)
  
# Save relative abundance data
  write.csv(y12, "OsmiaCC_Gen_bact_relabund.csv")
  
# Reorder x-axis  
  y12$combo_treat <- factor(y12$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot Genus by treatment
  OsmiaCC_gen_relabund_fungi <- ggplot(data = y12, aes(x = combo_treat, y = Abundance, fill = Genus)) + 
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
                                    ggtitle("B")
  OsmiaCC_gen_relabund_fungi
  
# Plot Genus for each sample
  ggplot(data = y12, aes(x = sampleID, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    facet_grid(~ combo_treat, 
               scale = "free", 
               space = "free") +
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Fungi")

## Differential abundance with raw data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# Remove patterns in tax_table   
  tax_table(ps4)[, colnames(tax_table(ps4))] <- gsub(tax_table(ps4)[, colnames(tax_table(ps4))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq::phyloseq_to_deseq2(ps4, ~ combo_treat)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means 
  geoMeans <- apply(counts(desq_obj), 1, gm_mean)
  
# Estimate size factors
  desq_dds <- estimateSizeFactors(desq_obj, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds <- DESeq(desq_dds, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# WN vs AN
  
# Extract results from differential abundance table for initial vs final provision
  WN_AN <- DESeq2::results(desq_dds, contrast = c("combo_treat", "WN", "AN"))
  
# Order differential abundances by their padj value
  WN_AN <- WN_AN[order(WN_AN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_AN_p05 <- WN_AN[(WN_AN$padj < alpha & !is.na(WN_AN$padj)), ]
  
# Check to see if any padj is below alpha
  WN_AN_p05
  
# AN vs CN

# Extract results from differential abundance table for initial vs final provision
  AN_CN <- DESeq2::results(desq_dds, contrast = c("combo_treat", "AN", "CN"))
  
# Order differential abundances by their padj value
  AN_CN <- AN_CN[order(AN_CN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AN_CN_p05 <- AN_CN[(AN_CN$padj < alpha & !is.na(AN_CN$padj)), ]
  
# Check to see if any padj is below alpha
  AN_CN_p05

# WN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  WN_CN <- DESeq2::results(desq_dds, contrast = c("combo_treat", "WN", "CN"))
  
# Order differential abundances by their padj value
  WN_CN <- WN_CN[order(WN_CN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_CN_p05 <- WN_CN[(WN_CN$padj < alpha & !is.na(WN_CN$padj)), ]
  
# Check to see if any padj is below alpha
  WN_CN_p05
  
## Differential abundance with rarefied data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# Convert from a phyloseq to a deseq obj
  desq_obj_rare <- phyloseq::phyloseq_to_deseq2(fung_rareps, ~ combo_treat)
  
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
