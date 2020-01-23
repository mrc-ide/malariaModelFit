# Prepare Griffin et al 2014 data

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

# Extract Incidence datasets
inc_data <- filter(data, !is.na(clin_k)) %>%
  mutate(numer = clin_k, denom = clin_p_years, case_detection = clinical_case_detection) %>%
  mutate(type = 1)
prev_data <- filter(data, !is.na(slide_k)) %>%
  mutate(numer = slide_k, denom = n, case_detection = clinical_case_detection) %>%
  mutate(type = 2)

study_data_2014 <- bind_rows(inc_data, prev_data) %>%
  filter(case_detection != -1) %>%
  rename(country_name = country, site_name = study) %>%
  mutate(reference = "TODO",
         study_index = study_index + 1,
         site_index = as.numeric(as.factor(paste(country_name, site_name, type)))) %>%
  select(country_name, site_name, reference, study_index, site_index, numer, denom, type,
         age0, age1, case_detection)

# Plot
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
pd
ggsave("data-raw/Griffin-et-al-2014-data/plot_data.png", pd, width = 11, height = 6)

# Add data to package
usethis::use_data(study_data_2014, overwrite = TRUE)
