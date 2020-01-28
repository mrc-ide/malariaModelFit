
# prepare_data_for_package.R
#
# Author: Pete Winskill
# Date: 2020-01-27
#
# Purpose:
# Read in original data used for Griffin et al. 2014 model fit. Process data
# into new format used within the malariaModelFit pacakge.
#
# ------------------------------------------------------------------

# load packages
library(dplyr)
library(ggplot2)

# Load, combine and name files
files <- list.files("data-raw/Griffin-et-al-2014-data/clin_data_used")
names <- unlist(strsplit(files, ".txt"))
data_files <- lapply(files, function(x){
  read.table(paste0("data-raw/Griffin-et-al-2014-data/clin_data_used/", x),
             sep = "\t")
})
names(data_files) <- names
data <- bind_rows(data_files, .id = "study")
colnames(data) <- c("study", "slide_k", "n", "pcr_k", "clin_k",
                    "clin_p_years", "sev_k", "sev_p_years", "age0", "age1",
                    "clin_group")

# Replace -1 missing/not relevant code with NA
data[data == -1] <- NA

# Link with study meta data
key <- read.csv("data-raw/Griffin-et-al-2014-data/study_info.csv", stringsAsFactors = FALSE)
data <- left_join(data, key, by = "study")

# Fix clinical_case_detection code
data$clinical_case_detection[data$clinical_case_detection == 6] <- 3

# Link with references
refs <- read.csv("data-raw/Griffin-et-al-2014-data/references.csv", stringsAsFactors = FALSE)
data <- left_join(data, refs, by = "study")

# Extract clinical incidence and prevalence datasets
inc_data <- filter(data, !is.na(clin_k)) %>%
  mutate(numer = clin_k, denom = clin_p_years, case_detection = clinical_case_detection) %>%
  mutate(type = 1)
prev_data <- filter(data, !is.na(slide_k)) %>%
  mutate(numer = slide_k, denom = n, case_detection = clinical_case_detection) %>%
  mutate(type = 2)

# Finalise data
study_data_2014 <- bind_rows(inc_data, prev_data) %>%
  filter(case_detection != -1) %>%
  rename(country_name = country, site_name = study) %>%
  mutate(study_index = match(study_index, unique(study_index)),
         site_index = as.numeric(as.factor(paste(type, site_name, country_name))),
         age_mid = age0 + (age1 - age0) / 2,
         age_bracket = case_when(age_mid < 2 ~ 1,
                                 age_mid >= 2 & age_mid < 5 ~ 2,
                                 age_mid >=5 & age_mid < 15 ~ 3,
                                 age_mid >= 15 ~ 4)) %>%
  select(country_name, site_name, reference, study_index, site_index, numer, denom, type,
         age0, age1, case_detection, age_bracket)

# Produce basic plot of all sites
pd <- ggplot(study_data_2014, aes(x = age0 + 0.5 * (age1 = age0),
                 y = numer / denom, shape = factor(type),
                 colour = country_name)) + 
  geom_line(alpha = 0.5, col = "black") +
  geom_point() + 
  xlab("Age") + 
  ylab("Y") + 
  ylim(0, NA) +
  scale_shape_discrete(labels = c("Inc", "Prev"), name = "Data type") + 
  facet_wrap(~ site_index, scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 6))

# Print plot to screen and to file
pd
ggsave("data-raw/Griffin-et-al-2014-data/plot_data.png", pd, width = 11, height = 6)

# Add data to package
usethis::use_data(study_data_2014, overwrite = TRUE)
