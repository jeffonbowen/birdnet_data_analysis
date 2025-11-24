

# ============================================================
# BirdNET-Pi QA + eBird BC barchart reference (NO API KEY)
# Uses your birds.db detections schema:
#   Date, Time, Sci_Name, Com_Name, Confidence, File_Name, ...
# ============================================================

library(DBI)
library(RSQLite)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(tidyr)

options(scipen = 999)

# ---------------------------
# USER SETTINGS (tweak later)
# ---------------------------

db_path <- "data/birds.db"
ebird_barchart_path <- "data/ebird_CA-BC__1900_2025_1_12_barchart.txt"

conf_low <- 0.70              # low-confidence flag threshold
conf_high_unlikely <- 0.90    # for rare/out-of-season, require this to pass
early_days_to_flag <- 2       # first N days of deployment flagged

season_present_thresh <- 0.001  # month considered "present" if freq >= this
# rarity thresholds based on annual mean relative frequency
thr_common   <- 0.02
thr_uncommon <- 0.005
thr_rare     <- 0.001
thr_veryrare <- 0.0002

# ---------------------------
# 1) LOAD BIRDNET DETECTIONS
# ---------------------------

con <- dbConnect(SQLite(), db_path)

det <- tbl(con, "detections") %>%
  select(Date, Time, Sci_Name, Com_Name, Confidence, File_Name) %>%
  collect()

dbDisconnect(con)

det <- det %>%
  mutate(
    Date = as.Date(Date),
    month_num = as.integer(substr(Date, 6, 7)),
    Com_Name = str_squish(Com_Name),
    Sci_Name = str_squish(Sci_Name)
  )

# Fix change in species names
det <- det %>%
  mutate(Com_Name = case_when(
    Com_Name == "Warbling Vireo" ~ "Western Warbling Vireo",
    Com_Name == "Herring Gull" ~ "American Herring Gull",
    Com_Name == "Pacific-slope Flycatcher" ~ "Western Flycatcher",
    Com_Name == "Yellow Warbler" ~ "Northern Yellow Warbler",
    TRUE ~ Com_Name
  ))

# helper: normalize names (remove apostrophes, trim whitespace)
normalize_name <- function(x) {
  x |>
    stringr::str_replace_all("['’]", "") |>  # remove straight & curly apostrophes
    stringr::str_squish()                    # trim + collapse whitespace
}

det <- det %>%
  mutate(Com_Name = normalize_name(Com_Name))

deploy_start <- min(det$Date, na.rm = TRUE)

# ----------------------------------------------------
# 2) LOAD + CLEAN EBIRD BC BARCHART REFERENCE
# ----------------------------------------------------
# File structure:
# - Metadata lines at top (Frequency of observations..., Sample Size...)
# - Then one row per taxon:
#     CommonName <tab> 48 numeric frequencies
# (You can see this in your file.) :contentReference[oaicite:1]{index=1}

raw_lines <- readLines(ebird_barchart_path)

# Find first species line by skipping until after "Sample Size:"
sample_idx <- which(str_detect(raw_lines, "^Sample Size:"))
if (length(sample_idx) == 0) stop("Couldn't find 'Sample Size:' line in barchart file.")
skip_n <- sample_idx[1] + 1

bc_ref <- read_tsv(
  ebird_barchart_path,
  skip = skip_n,
  col_names = FALSE,
  show_col_types = FALSE
)

# First column = common name, remaining columns = 48 period frequencies
names(bc_ref)[1] <- "Com_Name"
freq_cols <- names(bc_ref)[-1]

bc_ref <- bc_ref %>%
  mutate(
    Com_Name = str_squish(Com_Name)
  )

# Coerce frequencies to numeric (they come in as character sometimes)
bc_ref[freq_cols] <- lapply(bc_ref[freq_cols], as.numeric)

# ----------------------------------------------------
# 3) DERIVE MONTHLY + ANNUAL FREQUENCY FEATURES
# ----------------------------------------------------
# 48 bins = 4 per month (Jan..Dec). We'll average each 4-bin block.

# helper: indices for each month block
month_blocks <- split(seq_along(freq_cols), rep(1:12, each = 4))

# compute 12 month means
for (m in 1:12) {
  bc_ref[[paste0("m", m)]] <- rowMeans(bc_ref[, freq_cols[month_blocks[[m]]]], na.rm = TRUE)
}

month_mean_cols <- paste0("m", 1:12)

bc_ref <- bc_ref %>%
  mutate(
    annual_mean = rowMeans(across(all_of(freq_cols)), na.rm = TRUE),
    annual_max  = pmax(!!!syms(month_mean_cols), na.rm = TRUE),
    n_months_present = rowSums(across(all_of(month_mean_cols), ~ .x >= season_present_thresh),
                               na.rm = TRUE)
  )

# ----------------------------------------------------
# 4) ASSIGN RARITY CATEGORIES (BC-WIDE)
# ----------------------------------------------------
bc_ref <- bc_ref %>%
  mutate(
    rarity = case_when(
      annual_mean >= thr_common ~ "Common",
      annual_mean >= thr_uncommon ~ "Uncommon",
      annual_mean >= thr_rare ~ "Rare",
      annual_mean >= thr_veryrare ~ "Very rare",
      annual_mean > 0 ~ "Vagrant/Accidental",
      TRUE ~ "Listed but no freq"
    )
  )

# Month presence flags (TRUE/FALSE for each month)
for (m in 1:12) {
  bc_ref[[paste0("present_m", m)]] <- bc_ref[[paste0("m", m)]] >= season_present_thresh
}

present_cols <- paste0("present_m", 1:12)

# normalise names and add sort order
bc_ref <- bc_ref %>%
  mutate(Com_Name = normalize_name(Com_Name),
         sort_order = row_number())

# Have a look at the data
bc_ref_check <- bc_ref %>%
  select(Com_Name, sort_order, annual_mean, annual_max, n_months_present, rarity)
View(bc_ref_check)


# ----------------------------------------------------
# 5) JOIN DETECTIONS TO EBIRD REFERENCE
# ----------------------------------------------------
det_qa <- det %>%
  left_join(
    bc_ref %>% select(Com_Name, sort_order, annual_mean, annual_max, n_months_present, rarity, all_of(present_cols)),
    by = "Com_Name"
  ) %>%
  mutate(
    in_bc_list = !is.na(rarity)
  )

# ----------------------------------------------------
# 6) QA FLAGS (tiered)
# ----------------------------------------------------
det_qa <- det_qa %>%
  rowwise() %>%
  mutate(
    # A) Not in BC list at all (impossible)
    flag_not_in_bc = !in_bc_list,
    
    # B) Low confidence
    flag_low_conf = Confidence < conf_low,
    
    # C) Out-of-season using month presence
    flag_out_of_season = if_else(
      in_bc_list,
      !get(paste0("present_m", month_num)),
      TRUE  # if not in list, treat as out-of-season too
    ),
    
    # D) Rare/vagrant status
    flag_rare_status = in_bc_list & rarity %in% c("Rare", "Very rare", "Vagrant/Accidental", "Listed but no freq"),
    
    # E) Early deployment period
    flag_early_deploy = Date <= (deploy_start + days(early_days_to_flag - 1)),
    
    # Combined suspect rule
    flag_suspect =
      flag_not_in_bc |
      flag_early_deploy |
      (flag_out_of_season & Confidence < conf_high_unlikely) |
      (flag_rare_status & Confidence < conf_high_unlikely) |
      (flag_low_conf & flag_rare_status)
  ) %>%
  ungroup()

# ----------------------------------------------------
# 7) OUTPUTS FOR REVIEW + CLEAN DATASET
# ----------------------------------------------------

qa_flags_queue <- det_qa %>%
  filter(flag_suspect) %>%
  arrange(Date, Time, desc(Confidence))

qa_species_summary <- det_qa %>%
  group_by(Com_Name, Sci_Name, sort_order) %>%
  summarise(
    n = n(),
    max_conf = max(Confidence, na.rm = TRUE),
    mean_conf = mean(Confidence, na.rm = TRUE),
    rarity = first(rarity),
    annual_mean = first(annual_mean),
    n_months_present = first(n_months_present),
    prop_suspect = mean(flag_suspect),
    any_not_in_bc = any(flag_not_in_bc),
    any_out_of_season = any(flag_out_of_season),
    .groups = "drop"
  ) %>%
  arrange(desc(any_not_in_bc), desc(prop_suspect), desc(n))

det_clean <- det_qa %>%
  filter(!flag_suspect)



write_csv(det_qa, "data/detections_with_flags.csv")
write_csv(det_clean, "data/detections_clean.csv")
write_csv(qa_flags_queue, "data/qa_flags_queue.csv")
write_csv(qa_species_summary, "data/qa_species_summary.csv")
write_csv(bc_ref, "data/bc_ebird_reference_with_rarity.csv")

message("Saved outputs to /data:")
message(" - detections_with_flags.csv")
message(" - detections_clean.csv")
message(" - qa_flags_queue.csv")
message(" - qa_species_summary.csv")
message(" - bc_ebird_reference_with_rarity.csv")