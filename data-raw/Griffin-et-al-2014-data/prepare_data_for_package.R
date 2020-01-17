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
study_data_2014 <- left_join(data, key, by = "study")

# Visualise data to check
pd_prev <- study_data_2014 %>%
  filter(!is.na(slide_k)) %>%
  mutate(prev = slide_k / n,
         age = age0 + ((age1 - age0) / 2))
pd_clin <- study_data_2014 %>%
  filter(!is.na(clin_k)) %>%
  mutate(inc = clin_k / clin_p_years,
         age = age0 + ((age1 - age0) / 2))

p_prev <- ggplot(pd_prev, aes(x = age, y = prev, col = factor(study_index), group = 1)) +
  geom_line(alpha = 0.5, col = "black") +
  geom_point() + 
  facet_wrap(~ study) +
  theme_bw() +
  theme(legend.position="none",
        strip.text = element_text(size = 6))
p_inc <- ggplot(pd_clin, aes(x = age, y = inc, col = factor(study_index), group = 1)) +
  geom_line(alpha = 0.5, col = "black") +
  geom_point() + 
  facet_wrap(~ study) +
  theme_bw() +
  theme(legend.position="none",
        strip.text = element_text(size = 6))
ggsave("data-raw/Griffin-et-al-2014-data/plot_prev.png", p_prev, width = 10, height = 6)
ggsave("data-raw/Griffin-et-al-2014-data/plot_inc.png", p_inc, width = 8, height = 6)


# Add data to package
usethis::use_data(study_data_2014, overwrite = TRUE)
