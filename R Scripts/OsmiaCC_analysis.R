##### Project: Osmia Climate Change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose - Analyze microclimate, health, and life history data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(tidyverse) # Version 2.0.0
  library(lubridate) # Version 1.9.3
  library(lme4) # Version 1.1-35.1
  library(lmerTest) # Version 3.1-3
  library(pbkrtest) # Version 0.5.2
  library(emmeans) # Version 1.10.0
  library(effects) # Version 4.2-2
  library(ggplot2) # Version 3.4.3
  library(ggpubr) # Version 0.6.0
  library(cowplot) # Version 1.1.1
  library(knitr) # Version 1.45
  library(tibble) # Version 3.2.1
  library(survival) # Version 3.5-7
  library(ggsurvfit) # Version 1.0.0
  library(gtsummary) # Version 1.7.2
  library(coxme) # Version 2.2-18.1
  library(car) # Version 2.1-2
  library(survMisc) # Version 0.5.6
  library(multcomp) # Version 1.4-25
  library(multcompView) # Version 0.1-9

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
  all.sensors <- read.csv("all_sensors.csv")
  health <- read.csv("health - Data.csv")
  duration <- read.csv("life_history - Data.csv")
  mortality <- read.csv("mortality - Data.csv")

## Clean data ----

# Remove rows with incomplete data
  health <- na.omit(health)

# Subset to create dfs by sex
  males.health <- health[health$sex == "M", ]
  females.health <- health[health$sex == "F", ]
  males.duration <- duration[duration$sex == "M", ]
  females.duration <- duration[duration$sex == "F", ]
  males.mortality <- mortality[mortality$sex == "M", ]
  females.mortality <- mortality[mortality$sex == "F", ]

## Climate treatments ----
  
# Format date column
  all.sensors$Date <- as.Date(all.sensors$Date, format = "%Y-%m-%d")
  
# Remove dates after July 4, when the experiment stopped
  all.sensors <- all.sensors %>% filter(Date < '2023-07-04')

# Combine date and time columns
  all.sensors$DateTime <- as.POSIXct(paste(all.sensors$Date, all.sensors$Time),
                                     format = "%Y-%m-%d %I:%M:%S %p")  

# Manually order legend
  all.sensors$Sensor <- factor(all.sensors$Sensor, levels = c("Warm", "Ambient", "Cool"))

# Plot temperature & humidity
  treats <- ggplot(all.sensors, aes(x = DateTime, color = Sensor, group = Sensor)) +
                geom_line(aes(y = Temp, linetype = "Temperature")) + 
                geom_line(aes(y = Humidity, linetype = "Humidity")) +
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                theme(text = element_text(size = 16)) +
                xlab("Date") +
                scale_y_continuous(name = "Temperature (°C)", 
                                   sec.axis = sec_axis(trans = ~.*1, name = "Relative Humidity (%)")) +
                scale_color_manual(values = c("#C62828",  "#616161", "#1565C0"))
  treats
  
# Set Sensor reference group as ambient  
  all.sensors <- all.sensors %>% mutate(Sensor = relevel(Sensor, ref = "Ambient"))
  levels(all.sensors$Sensor)
  
# LMM of temperature using date as a random intercept to account for repeated measures
  gau.temp <- lmerTest::lmer(Temp ~ Sensor + (1|Date), data = all.sensors)
  summary(gau.temp)
  
# Pairwise comparisons with Tukey's HSD adjustment
  emmeans(gau.temp, pairwise ~ Sensor, adjust = "tukey")

# LMM of humidity using date as a random intercept to account for repeated measures
  gau.humidity <- lmerTest::lmer(Humidity ~ Sensor + (1|Date), data = all.sensors)
  summary(gau.humidity)
  
# Pairwise comparisons with Tukey's HSD adjustment
  emmeans(gau.humidity, pairwise ~ Sensor, adjust = "tukey")

# Metadata
  all.sensors %>%
    group_by(Sensor) %>%
    summarise(N = n(),
              Temp_mean = mean(Temp),
              Temp_se = sd(Temp)/sqrt(N),
              Temp_min = min(Temp),
              Temp_max = max(Temp),
              RH_mean = mean(Humidity),
              RH_se = sd(Humidity)/sqrt(N),
              RH_min = min(Humidity),
              RH_max = max(Humidity))
  
## Larval body mass ----

# Males  
  
# Subset df to include just treatments and response variable
  males.health.mass <- males.health %>%
    dplyr::select(bee, temp_treat, micro_treat, combo_treat, graft_stage, wet_mass_mg)

# Determine sample sizes, mean, and sd of males by treatment
  males.health.mass %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(wet_mass_mg), 
              SE = sd(wet_mass_mg)/sqrt(N))

# Change micro_treat and temp_treat to factors
  males.health.mass$micro_treat <- as.factor(males.health.mass$micro_treat)
  males.health.mass$temp_treat <- as.factor(males.health.mass$temp_treat)
  
# Set micro_treat reference group as sterile
  males.health.mass <- males.health.mass %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males.health.mass$micro_treat)

# Set temp_treat reference group as ambient  
  males.health.mass <- males.health.mass %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males.health.mass$temp_treat)
  
# LMM of larval biomass using grafting stage as a random intercept
  gau.mass.M <- lmerTest::lmer(wet_mass_mg ~ micro_treat * temp_treat + (1|graft_stage), data = males.health.mass)
 
# Check for normality with Q-Q plots and the Shapiro-Wilk test
  stats::qqnorm(resid(gau.mass.M, type = "pearson"))
  stats::qqline(resid(gau.mass.M, type = "pearson"))
  stats::shapiro.test(resid(gau.mass.M))
  
# LMM output
  summary(gau.mass.M)
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.mass.M, type = 2, ddf = "Kenward-Roger")
  
# Pairwise comparisons by temperature treatment with Tukey's HSD adjustment
  pairs(emmeans(gau.mass.M, "micro_treat", by = "temp_treat"), adjust = "tukey")
 
# Save p-values  
  stats.gau.mass.M <- tibble::tribble(
                                      ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                        "WS",     "WN",    0.02,      "*",
                                        "AS",     "WN",    0.04,      "*"
                                    )
  stats.gau.mass.M
   
# All pairwise comparisons with Tukey's HSD adjustment
  emmeans(gau.mass.M, pairwise ~ temp_treat * micro_treat, adjust = "tukey") 
  
# Reorder the x-axis
  males.health.mass$combo_treat <- factor(males.health.mass$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Italicize title
  bm.title <- expression(paste("(", italic("a"), ")"))
  
# Plot larval body mass by treatment
  bm.M <- ggplot(males.health.mass, aes(x = combo_treat, y = wet_mass_mg, color = combo_treat)) + 
            geom_boxplot(outlier.shape = NA,
                         width = 0.5,
                         position = position_dodge(width = 0.1)) + 
            geom_jitter(size = 1, 
                        alpha = 0.9) +
            theme_bw() +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = -0.17)) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            theme(text = element_text(size = 16)) +
            ylim(0, 20) +
            scale_color_manual(name = "Treatment", 
                               values = climate.colors) +
            labs(title = bm.title) +
            ylab("Larval body mass (mg)") +
            xlab("Treatment") +
            ggpubr::stat_pvalue_manual(stats.gau.mass.M,
                                       label = "p.adj.signif",
                                       y.position = 18,
                                       step.increase = 0.1,
                                       tip.length = 0.01)
  bm.M
  
# Females
  
# Subset df to include just treatments and response variable
  females.health.mass <- females.health %>%
    dplyr::select(bee, temp_treat, micro_treat, combo_treat, graft_stage, wet_mass_mg)
  
# Determine sample sizes, mean, and sd of males by treatment
  females.health.mass %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(wet_mass_mg), 
              SE = sd(wet_mass_mg)/sqrt(N))
  
# Change micro_treat and temp_treat to factors
  females.health.mass$micro_treat <- as.factor(females.health.mass$micro_treat)
  females.health.mass$temp_treat <- as.factor(females.health.mass$temp_treat)
  
# Set micro_treat reference group as sterile
  females.health.mass <- females.health.mass %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(females.health.mass$micro_treat)
  
# Set temp_treat reference group as ambient  
  females.health.mass <- females.health.mass %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(females.health.mass$temp_treat)
  
# LMM of larval biomass using grafting stage as a random intercept
  gau.mass.F <- lmerTest::lmer(wet_mass_mg ~ micro_treat * temp_treat + (1|graft_stage), data = females.health.mass)
  
# Check for normality with Q-Q plots and the Shapiro-Wilk test
  stats::qqnorm(resid(gau.mass.F, type = "pearson"))
  stats::qqline(resid(gau.mass.F, type = "pearson"))
  stats::shapiro.test(resid(gau.mass.F))
  
# LMM output
  summary(gau.mass.F)
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.mass.F, type = 2, ddf = "Kenward-Roger")
  
# Reorder the x-axis
  females.health.mass$combo_treat <- factor(females.health.mass$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot larval body mass by treatment
  bm.F <- ggplot(females.health.mass, aes(x = combo_treat, y = wet_mass_mg, color = combo_treat)) + 
              geom_boxplot(outlier.shape = NA,
                           width = 0.5,
                           position = position_dodge(width = 0.1)) + 
              geom_jitter(size = 1, 
                          alpha = 0.9) +
              theme_bw() +
              theme(legend.position = "none",
                    plot.title = element_text(hjust = -0.16)) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(text = element_text(size = 16)) +
              ylim(0, 20) +
              scale_color_manual(name = "Treatment", 
                                 values = climate.colors) +
              labs(title = bm.title) +
              ylab("Larval body mass (mg)") +
              xlab("Treatment")
  bm.F
  
## Proportion of larval body fat ----

# Males  
  
# Subset df to include just treatments and response variable
  males.health.fat <- males.health %>%
    dplyr::select(bee, temp_treat, micro_treat, combo_treat, graft_stage, prop_body_fat)
  
# Determine sample sizes, mean, and sd of males by treatment
  males.health.fat %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(prop_body_fat), 
              SD = sd(prop_body_fat))
  
# Change micro_treat and temp_treat to factors
  males.health.fat$micro_treat <- as.factor(males.health.fat$micro_treat)
  males.health.fat$temp_treat <- as.factor(males.health.fat$temp_treat)
  
# Set micro_treat reference group as sterile
  males.health.fat <- males.health.fat %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males.health.fat$micro_treat)
  
# Set temp_treat reference group as ambient  
  males.health.fat <- males.health.fat %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males.health.mass$temp_treat)
  
# LMM of larval fat content using grafting stage as a random intercept  
  gau.fat.M <- lmerTest::lmer(prop_body_fat ~ micro_treat * temp_treat + (1|graft_stage), data = males.health.fat)
  
# Check for normality with Q-Q plots and the Shapiro-Wilks test
  stats::qqnorm(resid(gau.fat.M, type = "pearson"))
  stats::qqline(resid(gau.fat.M, type = "pearson"))
  stats::shapiro.test(resid(gau.fat.M))
  
# LMM output
  summary(gau.fat.M)
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.fat.M, type = 2, ddf = "Kenward-Roger")
  
# Pairwise comparisons by temperature treatment with Tukey's HSD adjustment
  pairs(emmeans(gau.fat.M, "micro_treat", by = "temp_treat"), adjust = "tukey")
  
# All pairwise comparisons with Tukey's HSD adjustment
  emmeans(gau.fat.M, pairwise ~ temp_treat * micro_treat, adjust = "tukey")
  
# Save p-values  
  stats.gau.fat.M <- tibble::tribble(
                                     ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                       "WS",     "WN",    0.05,      "*",
                                       "AS",     "WN",    0.04,      "*"
                                 )
  stats.gau.fat.M
  
# Reorder the x-axis
  males.health.fat$combo_treat <- factor(males.health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Italicize title
  fat.title <- expression(paste("(", italic("b"), ")"))
  
# Plot the proportion of larval body fat by treatment
  fat.M <- ggplot(males.health.fat, aes(x = combo_treat, y = prop_body_fat, color = combo_treat)) + 
              geom_boxplot(outlier.shape = NA, 
                           width = 0.5, 
                           position = position_dodge(width = 0.1)) + 
              geom_jitter(size = 1, 
                          alpha = 0.9) +
              theme_bw() +
              theme(legend.position = "none",
                    plot.title = element_text(hjust = -0.15)) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(text = element_text(size = 16)) +
              ylim(0, 8) +
              scale_color_manual(name = "Treatment", 
                                 values = climate.colors) +
              labs(title = fat.title) +
              ylab("Proportion of body fat") + 
              xlab("Treatment") +
              ggpubr::stat_pvalue_manual(stats.gau.fat.M,
                                         label = "p.adj.signif",
                                         y.position = 6,
                                         step.increase = 0.1,
                                         tip.length = 0.01)
  fat.M
  
# Females  
  
# Subset df to include just treatments and response variable
  females.health.fat <- females.health %>%
    dplyr::select(bee, temp_treat, micro_treat, combo_treat, graft_stage, prop_body_fat)
  
# Determine sample sizes, mean, and sd of females by treatment
  females.health.fat %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(prop_body_fat), 
              SD = sd(prop_body_fat))
  
# Change micro_treat and temp_treat to factors
  females.health.fat$micro_treat <- as.factor(females.health.fat$micro_treat)
  females.health.fat$temp_treat <- as.factor(females.health.fat$temp_treat)
  
# Set micro_treat reference group as sterile
  females.health.fat <- females.health.fat %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males.health.fat$micro_treat)
  
# Set temp_treat reference group as ambient  
  females.health.fat <- females.health.fat %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(females.health.mass$temp_treat)
  
# LMM of larval fat content using grafting stage as a random intercept  
  gau.fat.F <- lmerTest::lmer(prop_body_fat ~ micro_treat * temp_treat + (1|graft_stage), data = females.health.fat)
  
# Check for normality with Q-Q plots and the Shapiro-Wilks test
  stats::qqnorm(resid(gau.fat.F, type = "pearson"))
  stats::qqline(resid(gau.fat.F, type = "pearson"))
  stats::shapiro.test(resid(gau.fat.F))
  
# LMM output
  summary(gau.fat.F)
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.fat.F, type = 2, ddf = "Kenward-Roger")
  
# Reorder the x-axis
  females.health.fat$combo_treat <- factor(females.health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot the proportion of larval body fat by treatment
  fat.F <- ggplot(females.health.fat, aes(x = combo_treat, y = prop_body_fat, color = combo_treat)) + 
              geom_boxplot(outlier.shape = NA, 
                           width = 0.5, 
                           position = position_dodge(width = 0.1)) + 
              geom_jitter(size = 1, 
                          alpha = 0.9) +
              theme_bw() +
              theme(legend.position = "none",
                    plot.title = element_text(hjust = -0.14)) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(text = element_text(size = 16)) +
              ylim(0, 8) +
              scale_color_manual(name = "Treatment", 
                                 values = climate.colors) +
              labs(title = fat.title) +
              ylab("Proportion of body fat") + 
              xlab("Treatment")
  fat.F
  
## Duration of developmental stages ----

# Males  
  
# Remove rows without complete data
  males.duration <- males.duration[!is.na(males.duration$days_instar2.5), ] 
  
# Determine sample sizes, mean, and sd of males by treatment
  males.duration %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(days_instar2.5), 
              SD = sd(days_instar2.5))
  
# Change micro_treat and temp_treat to factors
  males.duration$micro_treat <- as.factor(males.duration$micro_treat)
  males.duration$temp_treat <- as.factor(males.duration$temp_treat)
  
# Set micro_treat reference group as sterile
  males.duration <- males.duration %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males.duration$micro_treat)
  
# Set temp_treat reference group as ambient  
  males.duration <- males.duration %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males.duration$temp_treat)
  
# LMM of larval development using grafting stage as a random intercept 
  gau.dur.M <- lmerTest::lmer(days_instar2.5 ~ micro_treat * temp_treat + (1|graft_stage), data = males.duration)
  
# Check for normality with Q-Q plots and the Shapiro-Wilks test
  stats::qqnorm(resid(gau.dur.M, type = "pearson"))
  stats::qqline(resid(gau.dur.M, type = "pearson"))
  stats::shapiro.test(resid(gau.dur.M))
  
# GLMM of larval development using grafting stage as a random intercept and gamma distribution
  gam.dur.M <- lme4::glmer(days_instar2.5 ~ micro_treat * temp_treat + (1|graft_stage), family = Gamma, data = males.duration)

# Check for heteroscedasticity
  plot(gam.dur.M)
  
# GLMM output
  summary(gam.dur.M)
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.dur.M, type = 2, ddf = "Kenward-Roger")
  
# Pairwise comparisons by temperature treatment with Tukey's HSD adjustment
  pairs(emmeans(gam.dur.M, "micro_treat", by = "temp_treat"), adjust = "tukey")
  
#Pairwise comparisons by microbiome treatment with Tukey's HSD adjustment
  pairs(emmeans(gam.dur.M, "temp_treat", by = "micro_treat"), adjust = "tukey")
  
# All pairwise comparisons with Tukey's HSD adjustment
  emmeans(gam.dur.M, pairwise ~ temp_treat * micro_treat, adjust = "tukey")
  
# Save p-values
  stats.gam.dur.M <- tibble::tribble(
                                     ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                       "AS",     "WN",    0.03,      "*",
                                       "CS",     "WS",    0.01,      "*",
                                       "CS",     "WN",    0.0001,    "***",
                                       "WS",     "CN",    0.04,      "*",
                                       "CN",     "WN",    0.002,     "**"
                                   )
  stats.gam.dur.M
  
# Reorder the x-axis
  males.duration$combo_treat <- factor(males.duration$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN")) 
  
# Italicize title
  dur.title <- expression(paste("(", italic("c"), ")"))
  
# Plot larval duration by treatment
  dur.M <- ggplot(males.duration, aes(x = combo_treat, y = days_instar2.5, color = combo_treat)) + 
              geom_boxplot(outlier.shape = NA,
                           width = 0.5, 
                           position = position_dodge(width = 0.1)) + 
              geom_jitter(size = 1, 
                          alpha = 0.9) +
              theme_bw() +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(text = element_text(size = 16),
                    plot.title = element_text(hjust = -0.17)) +
              ylim(0, 25) +
              scale_color_manual(name = "Treatment", 
                                 values = climate.colors,
                                 labels = climate.labs) +             
              labs(title = dur.title) +
              ylab("Duration Larval instars II-V (days)") + 
              xlab("Treatment") +
              ggpubr::stat_pvalue_manual(stats.gam.dur.M,
                                         label = "p.adj.signif",
                                         y.position = 20,
                                         step.increase = 0.1,
                                         tip.length = 0.01)
  dur.M
  
# Females
  
# Remove rows without complete data
  females.duration <- females.duration[!is.na(females.duration$days_instar2.5), ] 
  
# Determine sample sizes, mean, and sd of males by treatment
  females.duration %>%
    group_by(combo_treat) %>%
    summarise(N = n(),
              Mean = mean(days_instar2.5), 
              SD = sd(days_instar2.5))
  
# Change micro_treat and temp_treat to factors
  females.duration$micro_treat <- as.factor(females.duration$micro_treat)
  females.duration$temp_treat <- as.factor(females.duration$temp_treat)
  
# Set micro_treat reference group as sterile
  females.duration <- females.duration %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(females.duration$micro_treat)
  
# Set temp_treat reference group as ambient  
  females.duration <- females.duration %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(females.duration$temp_treat)
  
# LMM of larval development using grafting stage as a random intercept 
  gau.dur.F <- lmerTest::lmer(days_instar2.5 ~ micro_treat * temp_treat + (1|graft_stage), data = females.duration)
  
# Check for normality with Q-Q plots and the Shapiro-Wilks test
  stats::qqnorm(resid(gau.dur.F, type = "pearson"))
  stats::qqline(resid(gau.dur.F, type = "pearson"))
  stats::shapiro.test(resid(gau.dur.F))
  
# ANOVA with Kenward-Roger approximation
  if(requireNamespace("pbkrtest", quietly = TRUE))
    anova(gau.dur.F, type = 2, ddf = "Kenward-Roger")
  
# Pairwise comparisons by temperature treatment with Tukey's HSD adjustment
  pairs(emmeans(gau.dur.F, "micro_treat", by = "temp_treat"), adjust = "tukey")
  
#Pairwise comparisons by microbiome treatment with Tukey's HSD adjustment
  pairs(emmeans(gau.dur.F, "temp_treat", by = "micro_treat"), adjust = "tukey")
  
# All pairwise comparisons with Tukey's HSD adjustment
  emmeans(gau.dur.F, pairwise ~ temp_treat * micro_treat, adjust = "tukey")
  
# Save p-values
  stats.gam.dur.F <- tibble::tribble(
                                     ~ group1, ~ group2, ~ p.adj, ~p.adj.signif,
                                       "AS",     "WS",     0.004,     "**",
                                       "AS",     "WN",     0.0006,    "**"
                                  )
  stats.gam.dur.F
  
# Reorder the x-axis
  females.duration$combo_treat <- factor(females.duration$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN")) 
  
# Plot larval duration by treatment
  dur.F <- ggplot(females.duration, aes(x = combo_treat, y = days_instar2.5, color = combo_treat)) + 
              geom_boxplot(outlier.shape = NA,
                           width = 0.5, 
                           position = position_dodge(width = 0.1)) + 
              geom_jitter(size = 1, 
                          alpha = 0.9) +
              theme_bw() +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(text = element_text(size = 16),
                    plot.title = element_text(hjust = -0.16)) +
              ylim(0, 25) +
              scale_color_manual(name = "Treatment", 
                                 values = climate.colors,
                                 labels = climate.labs) +             
              labs(title = dur.title) +
              ylab("Duration Larval instars II-V (days)") + 
              xlab("Treatment") +
              ggpubr::stat_pvalue_manual(stats.gam.dur.F,
                                         label = "p.adj.signif",
                                         y.position = 20,
                                         step.increase = 0.1,
                                         tip.length = 0.01)
  dur.F
  
## Mortality ----

# Males  
  
# How many bees from each treatment died within 48 h of grafting?
  males.mortality.graft <- males.mortality %>%
    filter(date > '6/8/2023') %>%
    group_by(combo_treat) %>%
    tally()
  males.mortality.graft
   
# Remove bees that died within 48 hr of grafting
  males.mortality48 <- males.mortality %>%
    filter(date < '6/8/2023')
  
# Determine sample sizes of males by treatment
  males.mortality.ss <- males.mortality48 %>%
    group_by(combo_treat) %>%
    tally()
  males.mortality.ss
  
# Add sample sizes per treatment (each originally had 30)
  males.mortality.ss$N <- c(29, 29, 29, 29, 30, 30)
  males.mortality.ss

# Add column and calculate percent mortality   
  males.mortality.ss$per_mort <- males.mortality.ss$n/males.mortality.ss$N
  males.mortality.ss

# Females
  
# How many bees from each treatment died within 48 h of grafting?
  females.mortality.graft <- females.mortality %>%
    filter(date > '6/9/2023') %>%
    group_by(combo_treat) %>%
    tally()
  females.mortality.graft
  
# Remove bees that died within 48 hr of grafting
  females.mortality48 <- females.mortality %>%
    filter(date < '6/9/2023')
  
# Determine sample sizes of females by treatment
  females.mortality.ss <- females.mortality48 %>%
    group_by(combo_treat) %>%
    tally()
  females.mortality.ss
  
# Add sample sizes per treatment (original ss: AN = 20, AS = 13, WN = 20, WS = 13)
  females.mortality.ss$N <- c(20, 13, 20, 13)
  females.mortality.ss
  
# Add column and calculate percent mortality   
  females.mortality.ss$per_mort <- females.mortality.ss$n/females.mortality.ss$N
  females.mortality.ss
  
## Survivorship analyses ----
# Resource: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# Resource: https://www.drizopoulos.com/courses/emc/basic_surivival_analysis_in_r#accelerated-failure-time-models
# NOTE: Status: 0 = survival to the fifth instar; 1 = death

# Both males and females  
  
# Recreate original df because NAs were removed above during analysis of larval development
  duration <- read.csv("life_history - Data.csv")
  
# Convert chr to date
  duration <- duration %>%
    mutate(
      date_nesting_start = ymd(date_nesting_start),
      date_nesting_end = ymd(date_nesting_end),
      date_graft = ymd(date_graft),
      date_instar1 = ymd(date_instar1),
      date_instar2 = ymd(date_instar2),
      date_instar5 = ymd(date_instar5),
      date_last_alive = ymd(date_last_alive)
    )
  
# Check to see if it worked  
  tibble(duration)
  
# Calculate the number of days bees survived
  duration <- duration %>%
    mutate(
      total_surv_days = as.duration(date_graft %--% date_last_alive / ddays(1))
    )

# Format chr to numeric
  duration$total_surv_days <- as.numeric(duration$total_surv_days)
  
# Check to see if it worked   
  head(duration)
  
# Create a survival object
  survival::Surv(duration$total_surv_days, duration$status)
  
# Fit the survival curve
  s1 <- ggsurvfit::survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat + sex, data = duration)
  summary(s1)
  
# Display mean survival time by combo_treat
  duration %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by temp_treat
  duration %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(N = n(),
              mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by micro_treat
  duration %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by micro_treat
  duration %>%
    filter(status == 1) %>%
    group_by(sex) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Log-rank test (when rho = 0) to compare survival times between groups to expected survival time
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat + sex, data = duration, rho = 0)
  
# Gehan-Wilcoxon test (when rho = 1) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat + sex, data = duration, rho = 1)
  
# Cox regression model to compare survival times between groups
  cox.mod1 <- coxme::coxme(Surv(total_surv_days, status) ~ sex + temp_treat * micro_treat + (1|graft_stage), data = duration)
  summary(cox.mod1)
  
# Test the proportional hazards assumption
  cz1 <- survival::cox.zph(cox.mod1)
  print(cz1)
  plot(cz1)

# Remove zeros in total_surv_days
  duration <- duration[duration$total_surv_days != 0, ]
  
# Accelerated failure time model
  aftm.all <- survival::survreg(Surv(total_surv_days, status) ~ sex + temp_treat + micro_treat + graft_stage, data = duration, dist = "weibull")
  summary(aftm.all)
  
# Use Kaplan-Meier estimator of the residuals to determine the fit of distribution
  
# Construct the residuals
  fitted.all <- aftm.all$linear.predictors
  res.all <- (log(aftm.all$y[, 1]) - fitted.all) / aftm.all$scale
  
# Fit regression model without covariates
  resKM.all <- survfit(Surv(res.all, status) ~ 1, data = duration)
  
# Plot residuals  
  plot(resKM.all, mark.time = FALSE, xlab = "AFTM Residuals", ylab = "Survival Probability")
  
# Superimpose the Weibull distribution
  xx <- seq(min(res.all), max(res.all), length.out = 35)
  yy <- exp(- exp(xx))
  lines(xx, yy, col = "red", lwd = 2)
  legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                         "Survival function of Extreme Value distribution"), 
         lty = c(1,2,1), col = c(1,1,2), bty = "n")
  
# Males only
  
# Recreate original df because NAs were removed above during analysis of larval development
  males.duration <- duration[duration$sex == "M", ]

# Convert chr to date
  males.duration <- males.duration %>%
    mutate(
      date_nesting_start = ymd(date_nesting_start),
      date_nesting_end = ymd(date_nesting_end),
      date_graft = ymd(date_graft),
      date_instar1 = ymd(date_instar1),
      date_instar2 = ymd(date_instar2),
      date_instar5 = ymd(date_instar5),
      date_last_alive = ymd(date_last_alive)
    )

# Check to see if it worked  
  tibble(males.duration)

# Calculate the number of days bees survived
  males.duration <- males.duration %>%
    mutate(
      total_surv_days = as.duration(date_graft %--% date_last_alive / ddays(1))
    )

# Format chr to numeric
  males.duration$total_surv_days <- as.numeric(males.duration$total_surv_days)

# Check to see if it worked   
  head(males.duration)

# NOTE: The analyses below includes bees that died within the first 48 h after grafting  
  
# Create a survival object
  survival::Surv(males.duration$total_surv_days, males.duration$status)

# Fit the survival curve
  s2 <- ggsurvfit::survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration)
  summary(s2)
  
# Display mean survival time by combo_treat
  males.duration %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by temp_treat
  males.duration %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))

# Display mean survival time by micro_treat
  males.duration %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Log-rank test (when rho = 0) to compare survival times between groups to expected survival time
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration, rho = 0)
  
# Gehan-Wilcoxon test (when rho = 1) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration, rho = 1)
  
# Change micro_treat and temp_treat to factors
  males.duration$micro_treat <- as.factor(males.duration$micro_treat)
  males.duration$temp_treat <- as.factor(males.duration$temp_treat)
  
# Set micro_treat reference group as sterile
  males.duration <- males.duration %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males.duration$micro_treat)
  
# Set temp_treat reference group as ambient  
  males.duration <- males.duration %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males.duration$temp_treat)
  
# Cox regression model to compare survival times between groups
  cox.mod2 <- coxme::coxme(Surv(total_surv_days, status) ~ temp_treat * micro_treat + (1|graft_stage), data = males.duration)
  summary(cox.mod2)
  
# Test the proportional hazards assumption
  cz2 <- survival::cox.zph(cox.mod2)
  print(cz2)
  plot(cz2)
  
# Accelerated failure time model
  aftm.M <- survival::survreg(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration, dist = "weibull")
  
# Construct the residuals
  fitted.M <- aftm.M$linear.predictors
  res.M <- (log(aftm.M$y[, 1]) - fitted.M) / aftm.M$scale
  
# Fit regression model without covariates
  resKM.M <- survfit(Surv(res.M, status) ~ 1, data = males.duration)
  
# Plot residuals  
  plot(resKM.M, mark.time = FALSE, xlab = "AFTM Residuals", ylab = "Survival Probability")
  
# Superimpose the Weibull distribution
  xx <- seq(min(res.M), max(res.M), length.out = 35)
  yy <- exp(- exp(xx))
  lines(xx, yy, col = "red", lwd = 2)
  legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                         "Survival function of Extreme Value distribution"), 
         lty = c(1,2,1), col = c(1,1,2), bty = "n")

# View model output  
  summary(aftm.M)
  
# Interpret coefficients 
  aftm.M.natural <- exp(1)^-0.0722
  aftm.M.natural <- 1 - aftm.M.natural
  aftm.M.natural
  
  aftm.M.cool <- exp(1)^-0.0649
  aftm.M.cool <- 1 - aftm.M.cool
  aftm.M.cool
  
  aftm.M.warm <- exp(1)^-0.2691
  aftm.M.warm <- 1 - aftm.M.warm
  aftm.M.warm
  
# Kaplan-Meier with all male bees
  OsmiaCC.KP.all.M <- ggsurvfit(s2) +
                        theme_classic() +
                        theme(legend.position = "right") +
                        theme(text = element_text(size = 16)) +
                        scale_color_manual(name = "Treatment", 
                                           values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                           labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                        scale_y_continuous(limits = c(0, 1)) +
                        scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                        labs(x = "Days",
                             y = "Survival probability")
  OsmiaCC.KP.all.M
  
# NOTE: The analyses below does NOT include bees that died within the first 48 h after grafting

# Remove bees that died within 48 h of grafting
  males.duration48 <- males.duration %>% 
    filter(date_last_alive > '2023-06-08')
  
# Change micro_treat and temp_treat to factors
  males.duration48$micro_treat <- as.factor(males.duration48$micro_treat)
  males.duration48$temp_treat <- as.factor(males.duration48$temp_treat)
  
# Set micro_treat reference group as sterile
  males.duration48 <- males.duration48 %>% mutate(micro_treat = relevel(micro_treat, ref = "natural"))
  levels(males.duration48$micro_treat)
  
# Set temp_treat reference group as ambient  
  males.duration48 <- males.duration48 %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males.duration48$temp_treat)
  
# Create a survival object
  survival::Surv(males.duration48$total_surv_days, males.duration48$status)
  
# Fit the survival curve
  s3 <- ggsurvfit::survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration48)
  summary(s3)
  
# Display mean survival time by combo_treat
  males.duration48 %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by temp_treat
  males.duration48 %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by micro_treat
  males.duration48 %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Log-rank test (when rho = 0) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration48, rho = 0)

# Gehan-Wilcoxon test (when rho = 1) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration48, rho = 1)
  
# Cox regression model to compare survival times between groups
  cox.mod3 <- coxme::coxme(Surv(total_surv_days, status) ~ temp_treat * micro_treat + (1|graft_stage), data = males.duration48)
  summary(cox.mod3)

# Test the proportional hazards assumption
  cz3 <- survival::cox.zph(cox.mod3)
  print(cz3)
  plot(cz3)
  
# Accelerated failure time model
  aftm.M.48 <- survival::survreg(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males.duration48, dist = "weibull")
  
# Construct the residuals
  fitted.M.48 <- aftm.M.48$linear.predictors
  res.M.48 <- (log(aftm.M.48$y[, 1]) - fitted.M.48) / aftm.M.48$scale
  
# Fit regression model without covariates
  resKM.M.48 <- survfit(Surv(res.M.48, status) ~ 1, data = males.duration48)
  
# Plot residuals  
  plot(resKM.M.48, mark.time = FALSE, xlab = "AFTM Residuals", ylab = "Survival Probability")
  
# Superimpose the Weibull distribution
  xx <- seq(min(res.M.48), max(res.M.48), length.out = 35)
  yy <- exp(- exp(xx))
  lines(xx, yy, col = "red", lwd = 2)
  legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                         "Survival function of Extreme Value distribution"), 
         lty = c(1,2,1), col = c(1,1,2), bty = "n")

# View model output
  summary(aftm.M.48)
  
# Interpret coefficients 
  aftm.M.48.sterile <- exp(1)^-0.0613
  aftm.M.48.sterile <- 1 - aftm.M.48.sterile
  aftm.M.48.sterile
  
  aftm.M.48.cool <- exp(1)^-0.1351
  aftm.M.48.cool <- 1 - aftm.M.48.cool
  aftm.M.48.cool
  
  aftm.M.48.warm <- exp(1)^-0.2099
  aftm.M.48.warm <- 1 - aftm.M.48.warm
  aftm.M.48.warm
  
# Kaplan-Meier without bees that died within 48 h of grafting
  OsmiaCC.KP.M.48 <- ggsurvfit(s3) +
                        theme_classic() +
                        theme(legend.position = "right") +
                        theme(text = element_text(size = 16)) +
                        scale_color_manual(name = "Treatment", 
                                           values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                           labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                        scale_y_continuous(limits = c(0, 1)) +
                        scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                        labs(x = "Days",
                             y = "Survival probability")
  OsmiaCC.KP.M.48
  
# Females only
  
# Recreate original df because NAs were removed above during analysis of larval development
  females.duration <- duration[duration$sex == "F", ]
  
# Convert chr to date
  females.duration <- females.duration %>%
    mutate(
      date_nesting_start = ymd(date_nesting_start),
      date_nesting_end = ymd(date_nesting_end),
      date_graft = ymd(date_graft),
      date_instar1 = ymd(date_instar1),
      date_instar2 = ymd(date_instar2),
      date_instar5 = ymd(date_instar5),
      date_last_alive = ymd(date_last_alive)
    )
  
# Check to see if it worked  
  tibble(females.duration)
  
# Calculate the number of days bees survived
  females.duration <- females.duration %>%
    mutate(
      total_surv_days = as.duration(date_graft %--% date_last_alive / ddays(1))
    )
  
# Format chr to numeric
  females.duration$total_surv_days <- as.numeric(females.duration$total_surv_days)
  
# Check to see if it worked   
  head(females.duration)
  
# NOTE: The analyses below includes bees that died within the first 48 h after grafting  
  
# Create a survival object
  survival::Surv(females.duration$total_surv_days, females.duration$status)
  
# Fit the survival curve
  s4 <- ggsurvfit::survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration)
  summary(s4)
  
# Display mean survival time by combo_treat
  females.duration %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by temp_treat
  females.duration %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by micro_treat
  females.duration %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Log-rank test (when rho = 0) to compare survival times between groups to expected survival time
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration, rho = 0)
  
# Gehan-Wilcoxon test (when rho = 1) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration, rho = 1)
  
# Cox regression model to compare survival times between groups
  cox.mod4 <- coxme::coxme(Surv(total_surv_days, status) ~ temp_treat * micro_treat + (1|graft_stage), data = females.duration)
  summary(cox.mod4)
  
# Test the proportional hazards assumption
  cz4 <- survival::cox.zph(cox.mod4)
  print(cz4)
  plot(cz4)
  
# Kaplan-Meier with all female bees
  OsmiaCC.KP.all.F <- ggsurvfit(s4) +
                        theme_classic() +
                        theme(legend.position = "right") +
                        theme(text = element_text(size = 16)) +
                        scale_color_manual(name = "Treatment", 
                                           values = c("#9E9E9E", "#616161", "#E57373", "#C62828"),
                                           labels = c('Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                        scale_y_continuous(limits = c(0, 1)) +
                        scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                        labs(x = "Days",
                             y = "Survival probability")
  OsmiaCC.KP.all.F
  
# NOTE: The analyses below does NOT include bees that died within the first 48 h after grafting
  
# Remove bees that died within 48 h of grafting
  females.duration48 <- females.duration %>% 
    filter(date_last_alive > '2023-06-09')
  
# Change micro_treat and temp_treat to factors
  females.duration48$micro_treat <- as.factor(females.duration48$micro_treat)
  females.duration48$temp_treat <- as.factor(females.duration48$temp_treat)
  
# Set micro_treat reference group as sterile
  females.duration48 <- females.duration48 %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(females.duration48$micro_treat)
  
# Set temp_treat reference group as ambient  
  females.duration48 <- females.duration48 %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(females.duration48$temp_treat)
  
# Create a survival object
  survival::Surv(females.duration48$total_surv_days, females.duration48$status)
  
# Fit the survival curve
  s5 <- ggsurvfit::survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration48)
  summary(s5)
  
# Display mean survival time by combo_treat
  females.duration48 %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by temp_treat
  females.duration48 %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Display mean survival time by micro_treat
  females.duration48 %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(N = n(),
              Mean = mean(total_surv_days),
              SE = sd(total_surv_days)/sqrt(N))
  
# Log-rank test (when rho = 0) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration48, rho = 0)
  
# Gehan-Wilcoxon test (when rho = 1) to compare survival times between groups
  survival::survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = females.duration48, rho = 1)
  
# Cox regression model to compare survival times between groups
  cox.mod5 <- coxme::coxme(Surv(total_surv_days, status) ~ temp_treat * micro_treat + (1|graft_stage), data = females.duration48)
  summary(cox.mod5)
  
# Test the proportional hazards assumption
  cz5 <- survival::cox.zph(cox.mod5)
  print(cz5)
  plot(cz5)
  
# Kaplan-Meier without bees that died within 48 h of grafting
  OsmiaCC.KP.F.48 <- ggsurvfit(s5) +
                        theme_classic() +
                        theme(legend.position = "right") +
                        theme(text = element_text(size = 16)) +
                        scale_color_manual(name = "Treatment", 
                                           values = c("#9E9E9E", "#616161", "#E57373", "#C62828"),
                                           labels = c('Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                        scale_y_continuous(limits = c(0, 1)) +
                        scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                        labs(x = "Days",
                             y = "Survival probability")
  OsmiaCC.KP.F.48
  
