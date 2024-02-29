##### Project: Osmia Climate Change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Create nice figures

### Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(patchwork) # Version 1.1.3

### Climate treatments ----

# Save plot
  ggsave("OsmiaCC_treatments.png", plot = treats, width = 10, height = 4, units = "in")

### Bee health and development ----

## Male bees

# Create plot
  OsmiaCC.fitness.M <- bm.M + fat.M + dur.M + plot_layout(ncol = 3)
  OsmiaCC.fitness.M

# Save plot  
  ggsave("OsmiaCC_fitness_M.png", plot = OsmiaCC.fitness.M, width = 15, height = 6, unit = "in")
  
## Female bees
  
# Create plot
  OsmiaCC.fitness.F <- bm.F + fat.F + dur.F + plot_layout(ncol = 3)
  OsmiaCC.fitness.F
  
# Save plot
  ggsave("OsmiaCC_fitness_F.png", plot = OsmiaCC.fitness.F, width = 15, height = 6, unit = "in")
  
### Kaplan-Meier ----
  
## All male bees
  
# Save plot
  ggsave("OsmiaCC_KP_all_M.png", plot = OsmiaCC.KP.all.M, width = 6, height = 4, unit = "in")

## Only male bees that did not die within 48 h
  
# Save plot
  ggsave("OsmiaCC_KP_M_48.png", plot = OsmiaCC.KP.M.48, width = 5, height = 3, unit = "in")
  
## All female bees
  
# Save plot
  ggsave("OsmiaCC_KP_all_F.png", plot = OsmiaCC.KP.all.F, width = 6, height = 4, unit = "in")

## Only female bees that did not die within 48 h
  
# Save plot
  ggsave("OsmiaCC_KP_F_48.png", plot = OsmiaCC.KP.F.48, width = 6, height = 4, unit = "in")
  
### Shannon Index ----

## Provisions without bees

# Create plot
  OsmiaCC.Shannon.NoBee <- OsmiaCC.Shannon.bact.NoBee + OsmiaCC.Shannon.fung.NoBee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Shannon.NoBee
  
# Save plot  
  ggsave("OsmiaCC_Shannon_NoBee.png", plot = OsmiaCC.Shannon.NoBee, width = 10, height = 4, units = "in")

## Provisions with bees
  
# Create plot
  OsmiaCC.Shannon.bee <- OsmiaCC.Shannon.bact.bee + OsmiaCC.Shannon.fung.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Shannon.bee
  
# Save plot  
  ggsave("OsmiaCC_Shannon_bee.png", plot = OsmiaCC.Shannon.bee, width = 10, height = 4, units = "in")
  
## Provisions with male bees
  
# Create plot
  OsmiaCC.Shannon.bee.M <- OsmiaCC.Shannon.bact.bee.M + OsmiaCC.Shannon.fung.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Shannon.bee.M
  
# Save plot  
  ggsave("OsmiaCC_Shannon_bee_M.png", plot = OsmiaCC.Shannon.bee.M, width = 10, height = 4, units = "in")
  
## Provisions with female bees
  
# Create plot
  OsmiaCC.Shannon.bee.F <- OsmiaCC.Shannon.bact.bee.F + OsmiaCC.Shannon.fung.bee.F + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Shannon.bee.F
  
# Save plot  
  ggsave("OsmiaCC_Shannon_bee_F.png", plot = OsmiaCC.Shannon.bee.F, width = 10, height = 4, units = "in")
  
### Simpson Index ----
  
## Provisions without bees

# Create plot
  OsmiaCC.Simpson.NoBee <- OsmiaCC.Simpson.bact.NoBee + OsmiaCC.Simpson.fung.NoBee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Simpson.NoBee
  
# Save plot  
  ggsave("OsmiaCC_Simpson_NoBee.png", plot = OsmiaCC.Simpson.NoBee, width = 10, height = 4, units = "in")
  
## Provisions with bees
  
# Create plot
  OsmiaCC.Simpson.bee <- OsmiaCC.Simpson.bact.bee + OsmiaCC.Simpson.fung.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Simpson.bee
  
# Save plot
  ggsave("OsmiaCC_Simpson_bee.png", plot = OsmiaCC.Simpson.bee, width = 10, height = 4, units = "in")
  
## Provisions with male bees
  
# Create plot
  OsmiaCC.Simpson.bee.M <- OsmiaCC.Simpson.bact.bee.M + OsmiaCC.Simpson.fung.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Simpson.bee.M
  
# Save plot
  ggsave("OsmiaCC_Simpson_bee_M.png", plot = OsmiaCC.Simpson.bee.M, width = 10, height = 4, units = "in")
  
## Provisions with female bees
  
# Create plot
  OsmiaCC.Simpson.bee.F <- OsmiaCC.Simpson.bact.bee.F + OsmiaCC.Simpson.fung.bee.F + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Simpson.bee.F
  
# Save plot
  ggsave("OsmiaCC_Simpson_bee_F.png", plot = OsmiaCC.Simpson.bee.F, width = 10, height = 4, units = "in")

### Observed Richness ----  
  
## Provisions without bees
  
# Create plot
  OsmiaCC.Observed.NoBee <- OsmiaCC.Observed.bact.NoBee + OsmiaCC.Observed.fung.NoBee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Observed.NoBee
  
# Save plot  
  ggsave("OsmiaCC_Observed_NoBee.png", plot = OsmiaCC.Observed.NoBee, width = 10, height = 4, units = "in")
  
## Provisions with bees
  
# Create plot
  OsmiaCC.Observed.bee <- OsmiaCC.Observed.bact.bee + OsmiaCC.Observed.fung.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Observed.bee
  
# Save plot
  ggsave("OsmiaCC_Observed_bee.png", plot = OsmiaCC.Observed.bee, width = 10, height = 4, units = "in")
  
## Provisions with male bees
  
# Create plot
  OsmiaCC.Observed.bee.M <- OsmiaCC.Observed.bact.bee.M + OsmiaCC.Observed.fung.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Observed.bee.M
  
# Save plot
  ggsave("OsmiaCC_Observed_bee_M.png", plot = OsmiaCC.Observed.bee.M, width = 10, height = 4, units = "in")

## Provisions with female bees
  
# Create plot
  OsmiaCC.Observed.bee.F <- OsmiaCC.Observed.bact.bee.F + OsmiaCC.Observed.fung.bee.F + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Observed.bee.F
  
# Save plot
  ggsave("OsmiaCC_Observed_bee_F.png", plot = OsmiaCC.Observed.bee.F, width = 10, height = 4, units = "in")

### Evenness ----  
  
## Provisions with bees
  
# Create plot
  OsmiaCC.Pielou.bee <- OsmiaCC.Pielou.bact.bee + OsmiaCC.Pielou.fung.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Pielou.bee
  
# Save plot
  ggsave("OsmiaCC_Pielou_bee.png", plot = OsmiaCC.Pielou.bee, width = 10, height = 4, units = "in")

## Male bees
  
# Create plot
  OsmiaCC.Pielou.M <- OsmiaCC.Pielou.bact.bee.M + OsmiaCC.Pielou.fung.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.Pielou.M
  
# Save plot
  ggsave("OsmiaCC_Pielou_bee_M.png", plot = OsmiaCC.Pielou.M, width = 10, height = 4, units = "in")
  
### PCoA with relative abundance data ----  

## Provisions with and without bees
  
# Create plot
  OsmiaCC.PCoA <- OsmiaCC.PCoA.bact + OsmiaCC.PCoA.fung + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA
  
# Save plot
  ggsave("OsmiaCC_PCoA.png", plot = OsmiaCC.PCoA, width = 10, height = 4, units = "in")
  
## Provisions without bees
  
# Create plot
  OsmiaCC.PCoA.NoBee <- OsmiaCC.PCoA.bact.NoBee + OsmiaCC.PCoA.fung.NoBee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.NoBee
  
# Save plot
  ggsave("OsmiaCC_PCoA_NoBee.png", plot = OsmiaCC.PCoA.NoBee, width = 10, height = 4, units = "in")
  
## Provisions bees
  
# Create plot
  OsmiaCC.PCoA.bee <- OsmiaCC.PCoA.bact.bee + OsmiaCC.PCoA.fung.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.bee
  
# Save plot
  ggsave("OsmiaCC_PCoA_bee.png", plot = OsmiaCC.PCoA.bee, width = 10, height = 4, units = "in")
  
## Provisions male bees
  
# Create plot
  OsmiaCC.PCoA.bee.M <- OsmiaCC.PCoA.bact.bee.M + OsmiaCC.PCoA.fung.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.bee.M
  
# Save plot
  ggsave("OsmiaCC_PCoA_bee_M.png", plot = OsmiaCC.PCoA.bee.M, width = 10, height = 4, units = "in")
  
## Provisions female bees
  
# Create plot
  OsmiaCC.PCoA.bee.F <- OsmiaCC.PCoA.bact.bee.F + OsmiaCC.PCoA.fung.bee.F + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.bee.F
  
# Save plot
  ggsave("OsmiaCC_PCoA_bee_F.png", plot = OsmiaCC.PCoA.bee.F, width = 10, height = 4, units = "in")
  
### Rarefaction curve ----
  
## Provisions with bees
  
# Create plot
  OsmiaCC.rare.bee <- OsmiaCC.rare.bact.bee + OsmiaCC.rare.fung.bee + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.rare.bee
  
# Save plot
  ggsave("OsmiaCC_rare_bee.png", plot = OsmiaCC.rare.bee, width = 5, height = 5, units = "in")
  
### PCoA with rarefied data  

## Provisions with bees  
  
# Create plot
  OsmiaCC.PCoA.rare.bee <- OsmiaCC.PCoA.bact.rare.bee + OsmiaCC.PCoA.fung.rare.bee + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.rare.bee

# Save plot
  ggsave("OsmiaCC_PCoA_rare_bee.png", plot = OsmiaCC.PCoA.rare.bee, width = 10, height = 4, units = "in")
  
## Provisions male bees
  
# Create plot
  OsmiaCC.PCoA.rare.bee.M <- OsmiaCC.PCoA.bact.rare.bee.M + OsmiaCC.PCoA.fung.rare.bee.M + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.rare.bee.M
  
# Save plot
  ggsave("OsmiaCC_PCoA_rare_bee_M.png", plot = OsmiaCC.PCoA.rare.bee.M, width = 10, height = 4, units = "in")
  
## Provisions female bees
  
# Create plot
  OsmiaCC.PCoA.rare.bee.F <- OsmiaCC.PCoA.bact.rare.bee.F + OsmiaCC.PCoA.fung.rare.bee.F + plot_layout(ncol = 2, nrow = 1)
  OsmiaCC.PCoA.rare.bee.F
  
# Save plot
  ggsave("OsmiaCC_PCoA_rare_bee_F.png", plot = OsmiaCC.PCoA.rare.bee.F, width = 10, height = 4, units = "in")
  
### Stacked community plots ----
  
## Control provisions
  
# Create plot
  OsmiaCC.stacked.controls <- OsmiaCC.gen.bact.controls + OsmiaCC.gen.fung.controls + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.controls
  
# Save plot
  ggsave("OsmiaCC_stacked_controls.png", plot = OsmiaCC.stacked.controls, width = 18, height = 12, units = "in")
  
## Provisions without bees
  
# Create plot
  OsmiaCC.stacked.NoBee <- OsmiaCC.gen.bact.NoBee + OsmiaCC.gen.fung.NoBee + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.NoBee
  
# Save plot
  ggsave("OsmiaCC_stacked_NoBee.png", plot = OsmiaCC.stacked.NoBee, width = 18, height = 12, units = "in")
  
## Provisions with bees
  
# Create plot
  OsmiaCC.stacked.bee <- OsmiaCC.gen.bact.bee + OsmiaCC.gen.fung.bee + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.bee
  
# Save plot
  ggsave("OsmiaCC_stacked_bee.png", plot = OsmiaCC.stacked.bee, width = 18, height = 12, units = "in")
  
## Provisions with male bees
  
# Create plot
  OsmiaCC.stacked.bee.M <- OsmiaCC.gen.bact.bee.M + OsmiaCC.gen.fung.bee.M + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.bee.M
  
# Save plot
  ggsave("OsmiaCC_stacked_bee_M.png", plot = OsmiaCC.stacked.bee.M, width = 18, height = 12, units = "in")
  
## Provisions with female bees
  
# Create plot
  OsmiaCC.stacked.bee.F <- OsmiaCC.gen.bact.bee.F + OsmiaCC.gen.fung.bee.F + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.bee.F
  
# Save plot
  ggsave("OsmiaCC_stacked_bee_F.png", plot = OsmiaCC.stacked.bee.F, width = 18, height = 12, units = "in")
  
## Top 15 genera in provisions with bees  
  
# Create plot
  OsmiaCC.stacked.15.bee <- OsmiaCC.15gen.relabund.bact + OsmiaCC.15gen.relabund.fung + plot_layout(ncol = 1, nrow = 2)
  OsmiaCC.stacked.15.bee
  
# Save plot
  ggsave("OsmiaCC_stacked_15_bee.png", plot = OsmiaCC.stacked.15.bee, width = 18, height = 12, units = "in")
  
### Arsenophonus ----
  
# Save plot
  ggsave("OsmiaCC_Arsenophonus_rel_abund.png", plot = arseno.rel.abund, width = 5, height = 3, units = "in")
  
