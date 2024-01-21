##### Project: Osmia climate change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose - Analyze health and life history data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(tidyverse) # Version 2.0.0
  library(lubridate) # Version 1.9.3
  library(ggplot2) # Version 3.4.3
  library(multcompView) # Version 0.1-9
  library(ggsignif) # Version 0.6.4
  library(patchwork) # Version 1.1.3
  library(cowplot) # Version 1.1.1
  library(knitr) # Version 1.45
  library(dplyr) # Version 1.1.3
  library(survival) # Version 3.5-7
  library(tibble) # Version 3.2.1
  library(ggsurvfit) # Version 1.0.0
  library(gtsummary) # Version 1.7.2
  library(car) # Version 2.1-2
  library(emmeans) # Version 1.8.9
  library(stats) # Version 4.3.1
  library(FSA) # Version 0.9.5
  library(lmtest) # Version 0.9-40
  library(survMisc) # Version 0.5.6

## Import data ----

# Import data
  all_sensors <- read.csv("all_sensors.csv")
  health <- read.csv("health - Data.csv")
  duration <- read.csv("life_history - Data.csv")
  mortality <- read.csv("mortality - Data.csv")

## Clean data ----

# Remove rows with no data
  health <- na.omit(health)

# Subset to create dfs by sex
  males_health <- health[health$sex == "M", ]
  males_duration <- duration[duration$sex == "M", ]
  males_mortality <- mortality[mortality$sex == "M", ]

## Analysis of temperature treatment data ----
  
# Format date column
  all_sensors$Date <- as.Date(all_sensors$Date, format = "%Y-%m-%d")
  
# Remove dates after July 4, when the experiment stopped
  all_sensors <- all_sensors %>% filter(Date < '2023-07-04')

# Combine date and time columns
  all_sensors$DateTime <- as.POSIXct(paste(all_sensors$Date, all_sensors$Time),
                                     format = "%Y-%m-%d %I:%M:%S %p")  

# Manually order legend
  all_sensors$Sensor <- factor(all_sensors$Sensor, levels = c("Warm", "Ambient", "Cool"))

# Plot temperature & humidity
  ggplot(all_sensors, aes(x = DateTime, color = Sensor, group = Sensor)) +
    geom_line(aes(y = Temp, linetype = "Temperature")) + 
    geom_line(aes(y = Humidity, linetype = "Humidity")) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(text = element_text(size = 16)) +
    xlab("Date") +
    scale_y_continuous(name = "Temperature (°C)", 
                       sec.axis = sec_axis(trans = ~.*1, name = "Relative Humidity (%)")) +
    scale_color_manual(values = c("#C62828",  "#616161", "#1565C0"))

## Larval body mass ----
  
# Subset df to include just treatments and response variable (wet body mass)
  males_health_mass <- males_health %>%
    select(bee, temp_treat, micro_treat, combo_treat, wet_mass_mg)
  
# Determine sample sizes by sex & treatment
  males_health_mass_ss <- males_health_mass %>%
    group_by(combo_treat) %>%
    tally()
  males_health_mass_ss
  
# ANOVA: Do any of the group means differ?
  males_mass_aov <- aov(wet_mass_mg ~ temp_treat * micro_treat, data = males_health_mass)

# Check for normality
  qqPlot(males_mass_aov$residuals, id = FALSE)  
  
# Shapiro-Wilk normality test
  shapiro.test(males_mass_aov$residuals)
  
# Levene's test to assess for equal variance
  leveneTest(males_health_mass$wet_mass_mg ~ males_health_mass$combo_treat)
  
# ANOVA output
  summary(males_mass_aov)
  
# Post-hoc: Which group means differ?  
  mass_tukey <- TukeyHSD(males_mass_aov)
  mass_tukey
  
# Linear contrasts
# Resource: https://bookdown.org/pingapang9/linear_models_bookdown/chap-contrasts.html#contrasts-and-dummy-coding

# What is the difference in wet body mass between sterile vs natural treatments? Dummy coding

# Ensure micro_treat and temp_treat are factors
  males_health_mass$micro_treat <- as.factor(males_health_mass$micro_treat)
  males_health_mass$temp_treat <- as.factor(males_health_mass$temp_treat)
  
# Set micro_treat reference group as sterile
  males_health_mass <- males_health_mass %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males_health_mass$micro_treat)

# Set temp_treat reference group as ambient  
  males_health_mass <- males_health_mass %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males_health_mass$temp_treat)

# Run linear model
  lm_mass1 <- lm(wet_mass_mg ~ micro_treat + temp_treat, data = males_health_mass)
  summary(lm_mass1)

# Display estimated marginal means for each pairwise comparison grouped by temperature treatment
  emmeans(lm_mass1, "micro_treat", by = "temp_treat")

# Display p-values for pairwise comparisons grouped by temperature treatment
  pairs(emmeans(lm_mass1, "micro_treat", by = "temp_treat"))

## Proportion of larval body fat ----

# Subset df to include just treatments and response variable (wet body mass)
  males_health_fat <- males_health %>%
    select(bee, temp_treat, micro_treat, combo_treat, prop_body_fat)
  
# Determine sample sizes by sex & treatment
  males_health_fat_ss <- males_health_fat %>%
    group_by(combo_treat) %>%
    tally()
  males_health_fat_ss
  
# ANOVA: Do any of the group means differ?
  males_fat_aov <- aov(prop_body_fat ~ temp_treat * micro_treat, data = males_health_fat)
  
# Check for normality
  qqPlot(males_fat_aov$residuals, id = FALSE)  
  
# Shapiro-Wilk normality test
  shapiro.test(males_fat_aov$residuals)
  
# Levene's test to assess for equal variance
  leveneTest(males_health_fat$prop_body_fat ~ males_health_fat$combo_treat)
  
# ANOVA output  
  summary(males_fat_aov)

# Post-hoc: Which group means differ?
  fat_tukey <- TukeyHSD(males_fat_aov)
  fat_tukey

# Linear contrasts
  
# What is the difference in total body fat between sterile vs natural treatments? Dummy coding
  
# Ensure micro_treat and temp_treat are factors
  males_health_fat$micro_treat <- as.factor(males_health_fat$micro_treat)
  males_health_fat$temp_treat <- as.factor(males_health_fat$temp_treat)
  
# Set micro_treat reference group as sterile
  males_health_fat <- males_health_fat %>% mutate(micro_treat = relevel(micro_treat, ref = "sterile"))
  levels(males_health_fat$micro_treat)
  
# Set temp_treat reference group as ambient  
  males_health_fat <- males_health_fat %>% mutate(temp_treat = relevel(temp_treat, ref = "ambient"))
  levels(males_health_mass$temp_treat)
  
# Run linear model
  lm_mass2 <- lm(prop_body_fat ~ micro_treat + temp_treat, data = males_health_fat)
  summary(lm_mass2)
  
# Display estimated marginal means for each pairwise comparison grouped by temperature treatment
  emmeans(lm_mass2, "micro_treat", by = "temp_treat")
  
# Display p-values for pairwise comparisons grouped by temperature treatment
  pairs(emmeans(lm_mass2, "micro_treat", by = "temp_treat"))

## Duration of developmental stages ----

# Determine sample sizes by sex & treatment
  males_duration_ss <- males_health %>%
    group_by(combo_treat) %>%
    tally()
  males_duration_ss
  
# ANOVA: Do any of the group means differ?
  males_dev_aov <- aov(days_instar2.5 ~ temp_treat + micro_treat, data = males_duration)
  
# Check for normality
  qqPlot(males_dev_aov$residuals, id = FALSE)  
  
# Shapiro-Wilk normality test
  shapiro.test(males_dev_aov$residuals)

# Levene's test to assess for equal variance
  leveneTest(males_duration$days_instar2.5 ~ males_duration$combo_treat)

# Kruskal-Wallis test: non-parametric one-way analysis of variance
  kruskal.test(days_instar2.5 ~ combo_treat, data = males_duration)

# Post-hoc: Dunn's test
  dunnTest(days_instar2.5 ~ combo_treat, data = males_duration, method = "bonferroni")

## Mortality ----

# Table of sample size per treatment
  males_mortality_ss <- males_mortality %>%
    group_by(combo_treat) %>%
    tally()
  males_mortality_ss

## Survivorship analysis ----
# Resource: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html  
# NOTE: Status: 0 = survival to the fifth instar; 1 = death
  
# Convert character to date
  males_duration <- males_duration %>%
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
  tibble(males_duration)

# Calculate the number of days bees survived
  males_duration <- males_duration %>%
    mutate(
      total_surv_days = as.duration(date_graft %--% date_last_alive / ddays(1))
    )

# Format chr to numeric
  males_duration$total_surv_days <- as.numeric(males_duration$total_surv_days)

# Check to see if it worked   
  head(males_duration)

# NOTE: The analyses below includes bees that died within the first 48 h after grafting  
  
# Create a survival object
  Surv(males_duration$total_surv_days, males_duration$status)

# Fit the survival curve
  s2 <- survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration)
  summary(s2)
  
# Display median survival time by combo_treat
  males_duration %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Display median survival time by temp_treat
  males_duration %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Display median survival time by micro_treat
  males_duration %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Log-rank test to compare survival times between groups (assumes risk of death to be same across time)
  survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration)

# Cox regression model to compare survival times between groups (allows risk of death to vary across time)
  cox_model1 <- coxph(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration)
  summary(cox_model1)
  
# Test the proportional hazards assumption
  par(mfrow = c(2, 1))

  cz1 <- cox.zph(cox_model1)
  print(cz1)
  plot (cz1)
  
# NOTE: The analyses below does NOT include bees that died within the first 48 h after grafting

# Remove bees that died within 48 h of grafting
  males_duration48 <- males_duration %>% 
    filter(date_last_alive > '2023-06-08')
  
# Create a survival object
  Surv(males_duration48$total_surv_days, males_duration48$status)
  
# Fit the survival curve
  s3 <- survfit2(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration48)
  summary(s3)
  
# Display median survival time by combo_treat
  males_duration48 %>%
    filter(status == 1) %>%
    group_by(combo_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Display median survival time by temp_treat
  males_duration48 %>%
    filter(status == 1) %>%
    group_by(temp_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Display median survival time by micro_treat
  males_duration48 %>%
    filter(status == 1) %>%
    group_by(micro_treat) %>%
    summarize(median_surv = median(total_surv_days))
  
# Log-rank test to compare survival times between groups (assumes risk of death to be same across time)
  survdiff(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration48)
  
# Cox regression model to compare survival times between groups (allows risk of death to vary across time)
  cox_model2 <- coxph(Surv(total_surv_days, status) ~ temp_treat + micro_treat, data = males_duration48)
  summary(cox_model2)

# Test the proportional hazards assumption
  par(mfrow = c(2, 1))
  
  cz2 <- cox.zph(cox_model2)
  print(cz2)
  plot (cz2)
  
## Plots ----  

# Body mass

# Reorder the x-axis
  males_health_mass$combo_treat <- factor(males_health_mass$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot larval body mass by treatment
  bm <- ggplot(males_health_mass, aes(x = combo_treat, y = wet_mass_mg,color = combo_treat)) + 
    geom_boxplot(outlier.shape = NA,
                 width = 0.5,
                 position = position_dodge(width = 0.1)) + 
    geom_jitter(size = 1, 
                alpha = 0.9) +
    theme_bw() +
    theme(legend.position = "none") + 
    theme(text = element_text(size = 16)) +
    ylim(0, 20) +
    scale_fill_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) +
    scale_x_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) + 
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
    ggtitle("A") +
    ylab("Larval body mass (mg)") + 
    xlab("Treatment")

# Total fat content
  
# Reorder the x-axis
  males_health_fat$combo_treat <- factor(males_health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# Plot the proportion of larval body fat by treatment
  fat <- ggplot(males_health_fat, aes(x = combo_treat, y = prop_body_fat, color = combo_treat)) + 
    geom_boxplot(outlier.shape = NA, 
                 width = 0.5, 
                 position = position_dodge(width = 0.1)) + 
    geom_jitter(size = 1, 
                alpha = 0.9) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(text = element_text(size = 16)) +
    ylim(0, 6) +
    scale_x_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) + 
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
    ggtitle("B") +
    ylab("Proportion of body fat") + 
    xlab("Treatment")
  
# Duration
  
# Reorder the x-axis
  males_health$combo_treat <- factor(males_health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN")) 
  
# Plot larval duration by treatment
  dur <- ggplot(males_duration, aes(x = combo_treat, y = days_instar2.5, color = combo_treat)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(size = 1, 
                alpha = 0.9) +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    ylim(5, 20) +
    scale_fill_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) +
    scale_x_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) + 
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +             
    ggtitle("C") +
    ylab("Duration Larval instars II-V (days)") + 
    xlab("Treatment")   
  
# Arrange health and life history plots
  OsmiaCC_fitness <- bm + fat + dur + plot_layout(ncol = 3)
  OsmiaCC_fitness
  
  ggsave("OsmiaCC_fitness.png", plot = OsmiaCC_fitness)

# Survivorship  
  
# Kaplan-Meier with all bees
  OsmiaCC_KP_all <- ggsurvfit(s2) +
                      theme_classic() +
                      theme(legend.position = "right") +
                      theme(text = element_text(size = 16)) +
                      scale_color_manual(name = "Treatment", 
                                         values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                         labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                      scale_y_continuous(limits = c(0, 1)) +
                      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                      labs(x = "Days", y = "Survival probability")
  OsmiaCC_KP_all
  
  ggsave("OsmiaCC_KP_all.png", plot = OsmiaCC_KP_all)
  
# Kaplan-Meier without bees that died within 48 h of grafting
  OsmiaCC_KP_48 <- ggsurvfit(s3) +
                      theme_classic() +
                      theme(legend.position = "right") +
                      theme(text = element_text(size = 16)) +
                      scale_color_manual(name = "Treatment", 
                                         values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                         labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                      scale_y_continuous(limits = c(0, 1)) +
                      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
                      labs(x = "Days", y = "Survival probability")
  OsmiaCC_KP_48
  
  ggsave("OsmiaCC_KP_48.png", plot = OsmiaCC_KP_48)
  
