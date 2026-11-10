# =============================================================================
#  Shoplifting in England: data assembly for the second ZISAR application.
#
#  Unit:     the 32,322 harmonised LSOA units of the solar application
#            (data/solar/solar_lsoa.rds, W_lsoa.rds and unit_lookup.rds,
#            written by solar_data.r), so that the two applications share
#            units, spatial weights and Census covariates
#  Outcome:  shoplifting offences recorded by the police from July 2016 to
#            June 2019 (36 months), from the police.uk street-level open data
#            (archive 2019-06.zip, 1.7 GB). Each offence is placed at an
#            anonymising snap point and carries the 2011 LSOA of that point;
#            offences with no LSOA (4%) are dropped
#  Covariates: employment by workplace from the Business Register and
#            Employment Survey 2018 (ONS, Nomis dataset NM_189_1) at LSOA
#            level, retail trade (SIC 2007 division 47) and all industries;
#            the BRES rounds small-area counts, so retail employment is also
#            banded. Census and IoD covariates come from solar_lsoa.rds
#  Police force: the force that recorded most of the unit's offences of all
#            types (British Transport Police excluded), for the force effects
#
#  Requires the duckdb command-line client and unzip.
#  Writes data/crime/shoplifting_lsoa.rds
# =============================================================================

library(Matrix)
library(dplyr)
library(readr)
library(tidyr)

raw_dir <- "../data/crime/raw"
out_dir <- "../data/crime"
dir.create(file.path(raw_dir, "bres"), recursive = TRUE, showWarnings = FALSE)

download_once <- function(url, dest) {
  if (!file.exists(dest)) download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

# -----------------------------------------------------------------------------
# 1. police.uk street-level files, Jul 2016 - Jun 2019: offences by month,
#    force, 2011 LSOA and crime type
# -----------------------------------------------------------------------------
zip_file <- download_once("https://policeuk-data.s3.amazonaws.com/archive/2019-06.zip",
                          file.path(raw_dir, "police_2019-06.zip"))
street_dir <- file.path(raw_dir, "street")
if (!dir.exists(street_dir))
  system2("unzip", c("-q", "-o", shQuote(zip_file), "'*-street.csv'", "-d", shQuote(street_dir)))
counts_csv <- file.path(raw_dir, "crime_counts_lsoa_month.csv")
if (!file.exists(counts_csv))
  system2("duckdb", c("-c", shQuote(sprintf(
    "SET threads = 4;
     COPY (SELECT Month AS month, \"Falls within\" AS force, \"LSOA code\" AS lsoa11,
                  \"Crime type\" AS crime_type, count(*) AS n
           FROM read_csv('%s/*/*-street.csv', header = true, all_varchar = true, union_by_name = true)
           GROUP BY ALL) TO '%s' (HEADER);", street_dir, counts_csv))))
cc <- read_csv(counts_csv, col_types = cols(n = "i", .default = "c"))
stopifnot(length(unique(cc$month)) == 36)
cat("offences:", sum(cc$n), " without LSOA:", sum(cc$n[is.na(cc$lsoa11)]),
    sprintf("(%.1f%%)", 100 * sum(cc$n[is.na(cc$lsoa11)]) / sum(cc$n)), "\n")
cc_all <- cc
cc <- filter(cc, startsWith(lsoa11, "E01"))           # England

# -----------------------------------------------------------------------------
# 2. BRES 2018: employment by workplace, retail (SIC 47) and all industries,
#    2011 LSOAs (Nomis API, 25,000 records per call)
# -----------------------------------------------------------------------------
bres_csv <- sapply(c(0, 25000, 50000), function(off) download_once(sprintf(paste0(
  "https://www.nomisweb.co.uk/api/v01/dataset/NM_189_1.data.csv?geography=TYPE298",
  "&industry=146800687,37748736&employment_status=4&measure=1&measures=20100&date=2018",
  "&select=geography_code,industry_code,obs_value,obs_status&recordoffset=%d"), off),
  file.path(raw_dir, "bres", sprintf("bres2018_%d.csv", off))))
bres <- bind_rows(lapply(bres_csv, read_csv, show_col_types = FALSE)) |>
  filter(startsWith(GEOGRAPHY_CODE, "E01")) |>
  transmute(lsoa11 = GEOGRAPHY_CODE, var = ifelse(INDUSTRY_CODE == "47", "retail_jobs", "jobs"),
            v = OBS_VALUE) |>
  pivot_wider(names_from = var, values_from = v)
stopifnot(nrow(bres) == 32844, !anyNA(bres$retail_jobs))

# -----------------------------------------------------------------------------
# 3. Units: the solar units, with the 2011 LSOA lookup of solar_data.r
# -----------------------------------------------------------------------------
solar <- readRDS("../data/solar/solar_lsoa.rds")
lk <- readRDS("../data/solar/unit_lookup.rds")           # LSOA11CD, unit
cc$unit <- lk$unit[match(cc$lsoa11, lk$LSOA11CD)]
bres$unit <- lk$unit[match(bres$lsoa11, lk$LSOA11CD)]
stopifnot(!anyNA(cc$unit), !anyNA(bres$unit))

shop <- cc |> filter(crime_type == "Shoplifting") |> count(unit, wt = n, name = "shoplifting")
offences <- cc |> count(unit, wt = n, name = "offences")
force <- cc |> filter(force != "British Transport Police") |>
  count(unit, force, wt = n) |> group_by(unit) |> slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |> select(unit, force)
bres_u <- bres |> group_by(unit) |> summarise(retail_jobs = sum(retail_jobs), jobs = sum(jobs), .groups = "drop")

crime <- solar |>
  select(unit, n_lsoa11, lad19cd, lad19nm, east, north, area_ha, households, pop, density,
         income_dep, sh_social, sh_flat, sh_private, sh_nssec12, sh_age65) |>
  left_join(shop, by = "unit") |> left_join(offences, by = "unit") |>
  left_join(force, by = "unit") |> left_join(bres_u, by = "unit") |>
  mutate(shoplifting = coalesce(shoplifting, 0L), offences = coalesce(offences, 0L),
         other_jobs = pmax(jobs - retail_jobs, 0),
         retail_band = cut(retail_jobs, c(-1, 0, 9, 24, 49, 99, 499, Inf),
                           labels = c("0", "1-9", "10-24", "25-49", "50-99", "100-499", "500+")),
         london = startsWith(lad19cd, "E09")) |>
  arrange(unit)
stopifnot(identical(crime$unit, solar$unit))
# units with no recorded offence take the modal force of their local authority;
# the City of London's five units are grouped with the Metropolitan Police
crime$force[crime$force == "City of London Police"] <- "Metropolitan Police Service"
lad_force <- tapply(crime$force, crime$lad19cd, function(x) names(which.max(table(x))))
crime$force <- coalesce(crime$force, unname(lad_force[crime$lad19cd]))
stopifnot(!anyNA(crime$force))
crime$force <- factor(crime$force)
saveRDS(crime, file.path(out_dir, "shoplifting_lsoa.rds"))

# -----------------------------------------------------------------------------
# 4. Descriptives
# -----------------------------------------------------------------------------
W <- readRDS("../data/solar/W_lsoa.rds")$queen
moran <- function(x) { z <- x - mean(x); as.numeric(crossprod(z, W %*% z) / crossprod(z)) }
y <- crime$shoplifting
cat("units:", nrow(crime), " shoplifting offences:", sum(y), " forces:", nlevels(crime$force), "\n")
cat(sprintf("mean %.2f variance %.1f median %d max %d; zeros %d (%.1f%%), outside London %.1f%%\n",
            mean(y), var(y), median(y), max(y), sum(y == 0), 100 * mean(y == 0), 100 * mean(y[!crime$london] == 0)))
print(quantile(y, c(.1, .25, .5, .75, .9, .95, .99)))
cat("\nzero share and counts by retail employment band:\n")
print(crime |> group_by(retail_band) |> summarise(units = n(), zero = round(mean(shoplifting == 0), 3),
      mean = round(mean(shoplifting), 1), median = median(shoplifting), .groups = "drop") |> as.data.frame())
cat(sprintf("\nMoran I: zero indicator %.3f, log(1 + y) %.3f, retail == 0 %.3f\n",
            moran(y == 0), moran(log1p(y)), moran(crime$retail_jobs == 0)))
cat("\nzero share by police force (lowest and highest):\n")
zf <- crime |> group_by(force) |> summarise(units = n(), zero = round(mean(shoplifting == 0), 3),
      per_1000 = round(1000 * sum(shoplifting) / sum(pop), 1), .groups = "drop") |> arrange(zero)
print(as.data.frame(rbind(head(zf, 5), tail(zf, 5))))
cat("\nmonths per force with fewer than half the force's median monthly offences:\n")
fm <- cc_all |> filter(!is.na(lsoa11), startsWith(lsoa11, "E01")) |> count(force, month, wt = n) |>
  group_by(force) |> filter(n < 0.5 * median(n), median(n) > 500) |> ungroup()
print(as.data.frame(fm))
