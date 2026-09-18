library(tidyverse)

#1. Read and prepare the text files containing CO2 release rates

files <- list.files(
  "Viktor_c_flux",
  pattern = "\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)

#Find header line

read_logger_file <- function(f) {
  lines <- readLines(f, warn = FALSE)
  
  # Find the real header
  header_line <- grep("DATE/TIME, DATA FORMAT", lines)
  
  # Keep everything from the header line onwards
  lines <- lines[header_line:length(lines)]
  
  # Remove wrapping quotes if any
  lines <- ifelse(
    str_starts(lines, '"') & str_ends(lines, '"'),
    str_sub(lines, 2, -2),
    lines
  )
  
  data <- read_delim(
    I(lines),
    delim = ",",
    trim_ws = TRUE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  # Add file name column
  data <- data %>% mutate(SourceFile = basename(f))
  
  return(data)
}


CO2_data <- map_dfr(files, read_logger_file) 

#Filter out all rows not containing CO2 data
CO2_data$CO2 <- as.numeric(CO2_data$CO2)
CO2_data <- CO2_data %>%
  filter(!is.na(CO2)) 

#Fix data structure
CO2_data <- CO2_data %>%
  mutate(
    TAIR = as.numeric(TAIR),
    AIRPRES = as.numeric(`AIR PRESSURE`),
    EventSec = period_to_seconds(hms(`EVENT TIME`)),
    Hour = str_extract(SourceFile, "(?<=_)[0-9]+(?=\\.txt)")
  ) %>%
  mutate(across(c(Hour, SourceFile, `EVENT DATE`), as.factor))

#Remove faulty measurements:
CO2_cleaned <- CO2_data %>%
  filter(!(`EVENT DATE` == "13/05/25" & SourceFile %in% c("bag1_9.txt",
                                                          "bag2_9.txt",
                                                          "bag3_13.txt",
                                                          "bag2_15.txt",
                                                          "bag1_16.txt",
                                                          "bag2_19.txt")))%>%
  filter(!(`EVENT DATE` == "25/05/25" & SourceFile %in% c("bag2_10.txt",
                                                          "bag3_11.txt",
                                                          "bag1_12.txt",
                                                          "bag2_12.txt",
                                                          "bag3_12.txt",
                                                          "bag2_13.txt",
                                                          "bag3_13.txt",
                                                          "bag3_14.txt",
                                                          "bag1_15.txt",
                                                          "bag2_15.txt",
                                                          "bag3_15.txt",
                                                          "bag1_16.txt",
                                                          "bag3_16.txt",
                                                          "bag2_17.txt",
                                                          "bag3_17.txt",
                                                          "bag1_18.txt",
                                                          "bag2_18.txt",
                                                          "bag3_18.txt",
                                                          "bag1_20.txt",
                                                          "bag1_22.txt",
                                                          "bag2_22.txt",
                                                          "bag3_22.txt")))%>%
  filter(!(`EVENT DATE` == "10/06/25" & SourceFile %in% c("bag2_10.txt",
                                                          "bag1_14.txt",
                                                          "bag2_15.txt",
                                                          "bag1_16.txt",
                                                          "bag3_18.txt")))




#Trim some measurements to start at the right time:
CO2_cleaned <- CO2_cleaned %>%
  filter(
    !(`EVENT DATE` == "10/06/25" & SourceFile == "bag1_12.txt") |
      EventSec >= 36690
  ) %>%
  filter(
    !(`EVENT DATE` == "13/05/25" & SourceFile == "bag2_16.txt") |
      EventSec >= 50306
  ) %>%
  filter(
    !(`EVENT DATE` == "25/05/25" & SourceFile == "bag3_20.txt") |
      EventSec >= 65275
  ) %>%
  filter(
    !(`EVENT DATE` == "25/05/25" & SourceFile == "bag3_11.txt") |
      EventSec >= 32030
  )

#2. Calculate ml CO2 per minute between 1000-2000ppm

#Get time in seconds
CO2_filtered <- CO2_cleaned %>%
  filter(CO2 >= 1000, CO2 <= 2000)

#Calculate release rate in ppm per minute
CO2_min <- CO2_filtered %>%
  group_by(`EVENT DATE`, SourceFile) %>%
  summarize(
    TotalRelease = max(CO2) - min(CO2), #How much was released
    TimeSec = max(EventSec) - min(EventSec), #During how many seconds
    ReleasePerMin = TotalRelease/TimeSec * 60,
    Temp = mean(TAIR),
    Pressure = mean(AIRPRES)
  )

#Extract hours
CO2_min <- CO2_min %>%
  mutate(
    Hour = str_extract(SourceFile, "(?<=_)[0-9]+(?=\\.txt)")
  ) %>%
  mutate(Hour = as.numeric(Hour))  

#Calculate ppm to ml
#CO2ml = Volume of chamber in mm * CO2/10^6 * Pressure * 273.15/Temp in Kelvin
#OBS: Pressure and temperature already accounted for in raw CO2 measurements!
#Thus CO2ppm = CO2 * Pressure * 273.15/Temp in Kelvin
#CO2ml = Volume of chamber in mm * CO2ppm / 10^6
  
#Calculate box volume in cubic meters 33x33x34cm =>
volume <- 0.33 * 0.33 * 0.34 

#Convert to milliliters (1 kubic meter = 1000 liters = 10^6 milliliters)
volume_ml <- volume * 10^6

#Convert to ml per minute
CO2_min <- CO2_min %>%
  mutate(
    mlpermin = volume_ml * ReleasePerMin/10^6
  ) 

write.csv(CO2_min, "CO2_per_minute.csv")