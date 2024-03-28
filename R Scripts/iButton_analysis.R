##### Project: Osmia climate change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Use temp. data from iButton loggers (Little Bear Creek, Logan, UT 2022-2023) to create a daily temperature regime

# Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(tidyverse) # Version 2.0.0
  library(lubridate) # Version 1.9.3
  library(dplyr) # Version 1.1.4
  library(readr) # Version 2.1.5
  library(xts) # Version 0.13.2
  library(ggplot2) # Version 3.4.4

# Import and format data

# Import data
  all_loggers <- list.files(path = "iButton_loggers", full.names = TRUE) %>% 
    lapply(read_csv) %>% 
      bind_rows
  
# Rename column headers
  colnames(all_loggers) <- c("Date", "Time", "Unit", "Temp")
  
# Format Date column
  all_loggers$Date <- as.Date(all_loggers$Date, format = "%m/%d/%y")
  
# Remove days when the iButton loggers were in the lab
  all_loggers <- all_loggers %>% filter(Date > "2022-06-13" & Date < "2022-10-29")

# Analyze data by day, week & month ----
  
# Group data by day, and then calculate the mean, maximum and minimum temperature
  daily_temps <- all_loggers %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(Mean = mean(Temp),
                     Max = max(Temp),
                     Min = min(Temp))
  daily_temps
  
# Add week in a new column
  daily_temps$Week <- strftime(daily_temps$Date, format = "%V")  
  
# Group data by week, and then calculate the mean, maximum and minimum temperature
  weekly_temps <- daily_temps %>%
    dplyr::group_by(Week) %>%
    dplyr::summarise(Mean = mean(Mean),
                     Min = mean(Min),
                     Max = mean(Max))
  weekly_temps
  
# Add month in a new column
  daily_temps$Month <- strftime(daily_temps$Date, format = "%m")
  
# Group data by month, and then calculate the mean, maximum and minimum temperature
  monthly_temps <- daily_temps %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(Mean = mean(Mean),
                     Min = mean(Min),
                     Max = mean(Max))
  monthly_temps

# Hourly mean temps by day ----
  
# Format into a xts and add column combining Date and Hour
  all_loggers.xts <- xts(all_loggers$Temp, as.POSIXct(paste(all_loggers$Date, all_loggers$Time), format= "%Y-%m-%d %H:%M:%S"))
  
# Check to see if it worked
  head(all_loggers.xts)
  
# Using data from all loggers, calculate the means for each hour
  means <- xts::period.apply(all_loggers.xts, endpoints(all_loggers.xts, "hours"), mean)

# Check to see if it worked  
  head(means)
  
# Align time stamps to the beginning of the hour
  align.time.down <- function(x,n){index(x) = index(x)-n; align.time(x,n)}
  means.rounded <- align.time.down(means, 60*60)
  
# Check to see if it worked
  head(means.rounded)
  
# Plot data
  plot(means.rounded, las = 1, main = "Hourly Avg Temperatures")

# Format into a data frame
  means.rounded.new <- data.frame(Date = index(means.rounded), coredata(means.rounded))
  
# Extract date; put into a new column
  means.rounded.new$Day <- as.Date(means.rounded.new$Date)
  
# Extract time; put into a new column
  means.rounded.new$Time <- format(as.POSIXct(means.rounded.new$Date), format = "%H:%M%S")
  
# Add month in a new column
  means.rounded.new$Month <- strftime(means.rounded.new$Day, format = "%m")
  
# Change the column name
  colnames(means.rounded.new)[2] = "Temp"  

# Hourly mean temperatures by month ----
  
# Filter data by month
  june <- means.rounded.new %>%
    filter(Month == "06")
  
# Group data by hour, and then calculate the mean, maximum and minimum temperature
  june.temps <- june %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(Mean = mean(Temp),
              Max = max(Temp),
              Min = min(Temp))
  june.temps
  
# Save data
  write.csv(june.temps, "june_temps.csv", row.names = FALSE)
  
# Plot mean temp by hour
    ggplot(june.temps, aes(x = Time, y = Mean)) +
        geom_point() + 
        ggtitle("June") +
        ylab("Temperature (°C)") +
        scale_x_discrete(labels = 0:23) +
        ylim(-5, 35) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank())
    
# Filter data by month
  july <- means.rounded.new %>%
    filter(Month == "07")
  
# Group data by hour, and then calculate the mean, maximum and minimum temperature
  july.temps <- july %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(Mean = mean(Temp),
                     Max = max(Temp),
                     Min = min(Temp))
  july.temps
  
# Plot mean temp by hour
    ggplot(july.temps, aes(x = Time, y = Mean)) +
        geom_point() + 
        ggtitle("July") +
        ylab("Temperature (°C)") +
        scale_x_discrete(labels = 0:23) +
        ylim(-5, 35) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank())
  
# Filter data by month
  august <- means.rounded.new %>%
    filter(Month == "08")
  
# Group data by hour, and then calculate the mean, maximum and minimum temperature
  august.temps <- august %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(Mean = mean(Temp),
                     Max = max(Temp),
                     Min = min(Temp))
  
# Plot mean temp by hour
    ggplot(august.temps, aes(x = Time, y = Mean)) +
        geom_point() + 
        ggtitle("August") +
        ylab("Temperature (°C)") +
        scale_x_discrete(labels = 0:23) +
        ylim(-5, 35) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank())    
  
# Filter data by month
  september <- means.rounded.new %>%
    filter(Month == "09")
  
# Group data by hour, and then calculate the mean, maximum and minimum temperature
  september.temps <- september %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(Mean = mean(Temp),
                     Max = max(Temp),
                     Min = min(Temp))  
  
# Plot mean temp by hour
    ggplot(september.temps, aes(x = Time, y = Mean)) +
        geom_point() + 
        ggtitle("September") +
        ylab("Temperature (°C)") +
        scale_x_discrete(labels = 0:23) +
        ylim(-5, 35) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank())   
  
# Filter data by month
  october <- means.rounded.new %>%
    filter(Month == "10")
  
# Group data by hour, and then calculate the mean, maximum and minimum temperature
  october.temps <- october %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(Mean = mean(Temp),
                     Max = max(Temp),
                     Min = min(Temp))    

# Plot mean temp by hour
    ggplot(october.temps, aes(x = Time, y = Mean)) +
        geom_point() + 
        ggtitle("October") +
        ylab("Temperature (°C)") +
        scale_x_discrete(labels = 0:23) +
        ylim(-5, 35) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(panel.grid.minor = element_blank())
