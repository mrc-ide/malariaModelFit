# Prepare Griffin et al 2014 data

library(dplyr)

# Load, combine and name files
files <- list.files("data-raw/Griffin-et-al-2014-data/clin_data_used")
names <- unlist(strsplit(files, ".txt"))
data_files <- lapply(files, function(x){
  read.table(paste0("data-raw/Griffin-et-al-2014-data/clin_data_used/", x),
             sep = "\t")
})
names(data_files) <- names
data <- bind_rows(data_files, .id = "study")
colnames(data) <- c("study", "slide_k", "slide_n", "pcr_k", "clin_k",
                    "clin_p_years", "sev_k", "sev_p_years", "age0", "age1",
                    "clin_group")

# Replace -1 missing/not relevant code with NA
data[data == -1] <- NA

# Link with study meta data
key <- read.csv("data-raw/Griffin-et-al-2014-data/study_info.csv", stringsAsFactors = FALSE)
study_data_2014 <- left_join(data, key, by = "study")

# Add data to package
usethis::use_data(study_data_2014)
