##### Project: Osmia climate change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose - Analyze health and life history data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  library(multcompView)
  library(ggsignif)
  library(gridExtra)
  library(cowplot)
  library(knitr)
  library(dplyr)
  library(survival)
  library(tibble)
  library(ggsurvfit)
  library(gtsummary)

## Import & format data ----

# Import data
  all_sensors <- read.csv("all_sensors.csv")
  health <- read.csv("health - Raw Data.csv")
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

# Remove dates after July 4, when the experiment stopped
  all_sensors <- all_sensors %>% filter(date < '2023-07-04')

# Format date column
  all_sensors$date <- as.Date(all_sensors$date, format = "%Y-%m-%d")

# Combine date and time columns
  all_sensors$DateTime <- as.POSIXct(paste(all_sensors$date, all_sensors$time),
                                     format = "%Y-%m-%d %I:%M:%S %p")  

# Manually order legend
  all_sensors$sensor <- factor(all_sensors$sensor, levels = c("warm", "ambient", "cool"))

# Plot temperature & humidity
  ggplot(all_sensors, aes(x = DateTime, color = sensor, group = sensor)) +
    geom_line(aes(y = temp, 
                  linetype = "Temperature")) + 
    geom_line(aes(y = humidity, 
                  linetype = "Humidity")) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(text = element_text(size = 16)) +
    xlab("Date") +
    scale_y_continuous(name = "Temperature (°C)", 
                       sec.axis = sec_axis(trans = ~.*1, name = "Relative Humidity (%)")) +
    scale_color_manual(values = c("#C62828",  "#616161", "#1565C0"))  

## Larval body mass ----

# Determine sample sizes by sex & treatment
  males_health_ss <- males_health %>%
    group_by(combo_treat) %>%
    tally()
  males_health_ss

# Subset df to include just treatments and response variable (wet body mass)
  males_health_mass <- males_health %>%
    select(bee, temp_treat, micro_treat, combo_treat, wet_mass_mg)

# ANOVA: Do any of the group means differ?  
  males_mass_lm <- lm(wet_mass_mg ~ temp_treat * micro_treat, data = males_health_mass)
  males_mass_aov <- aov(males_mass_lm)
  summary(males_mass_aov)  

# Post-hoc: Which group means differ?  
  mass_tukey <- TukeyHSD(males_mass_aov)
  mass_tukey

# Reorder the x-axis
  males_health_mass$combo_treat <- factor(males_health_mass$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))  

# Plot larval body mass by treatment
  ggplot(males_health_mass, aes(x = combo_treat, y = wet_mass_mg,color = combo_treat)) + 
    geom_boxplot(outlier.shape = NA, 
                 width = 0.5, 
                 position = position_dodge(width = 0.1)) + 
    geom_jitter(size = 1, 
                alpha = 0.9) +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    ylim(0, 20) +
    scale_fill_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) +
    scale_x_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) + 
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
    ggtitle("") +
    ylab("Larval body mass (mg)") + 
    xlab("Treatments")  

## Proportion of larval body fat ----

# Subset df to include just treatments and response variable (wet body mass)
  males_health_fat <- males_health %>%
    select(bee, temp_treat, micro_treat, combo_treat, prop_body_fat)

# ANOVA: Do any of the group means differ?
  males_fat_lm <- lm(prop_body_fat ~ temp_treat * micro_treat, data = males_health_fat)
  males_fat_aov <- aov(males_fat_lm)
  summary(males_fat_aov)

# Post-hoc: Which group means differ?
  fat_tukey <- TukeyHSD(males_fat_aov)
  fat_tukey

# Reorder the x-axis
  males_health_fat$combo_treat <- factor(males_health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# Plot the proportion of larval body fat by treatment
  ggplot(males_health_fat, aes(x = combo_treat, y = prop_body_fat, color = combo_treat)) + 
    geom_boxplot(outlier.shape = NA, 
                 width = 0.5, 
                 position = position_dodge(width = 0.1)) + 
    geom_jitter(size = 1, 
                alpha = 0.9) +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    ylim(0, 6) +
    scale_x_discrete(labels = c("CS", "CN", "AS", "AN", "WS", "WN")) + 
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
    ggtitle("") +
    ylab("Proportion of body fat") + 
    xlab("Treatment")

## Duration of developmental stages ----

# ANOVA: Do any of the group means differ?
  males_dev_lm <- lm(days_instar2.5 ~ temp_treat*micro_treat, data = males_duration)
  anova_dev <- aov(males_dev_lm)
  summary(anova_dev)

# Post-hoc: Which group means differ?
  ph_dev_tukey <- TukeyHSD(anova_dev)
  ph_dev_tukey

# Reorder the x-axis
  males_health$combo_treat <- factor(males_health$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN")) 

# Plot larval duration by treatment
  ggplot(males_duration, aes(x = combo_treat, y = days_instar2.5, color = combo_treat)) + 
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
    ggtitle("") +
    ylab("Duration Larval instars II-V (days)") + 
    xlab("Treatment") 

## Mortality

# Table of sample size per treatment
  males_mortality_ss <- males_mortality %>%
    group_by(combo_treat) %>%
    tally()

# View table  
  males_mortality_ss

## Survivorship analysis ----
# Resource: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html  
# NOTE: Status: 0 = survival to the fifth instar; 1 = death

# Convert chr to date
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

# Create a survival object
  Surv(males_duration$total_surv_days, males_duration$status)

# Create a survival curve
  survfit(Surv(total_surv_days, status) ~ combo_treat, data = males_duration)

# Fit the line
  s2 <- survfit2(Surv(total_surv_days, status) ~ combo_treat, data = males_duration)

# Plot Kaplan-Meier
  ggsurvfit(s2) +
    theme_classic() +
    theme(legend.position = "right") +
    theme(text = element_text(size = 16)) +
    scale_color_manual(name = "Treatment", 
                       values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                       labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
    labs(x = "Days", y = "Survival probability")

# Cox regression model
  coxph(Surv(total_surv_days, status) ~ combo_treat, data = males_duration) 

# Compare survival times between treatments
  survdiff(Surv(total_surv_days, status) ~ combo_treat, data = males_duration)
