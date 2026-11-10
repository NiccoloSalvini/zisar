# =============================================================================
#  Domestic solar PV adoption in England: data assembly for the ZISAR
#  application.
#
#  Unit:     Lower Layer Super Output Area (LSOA), England. Ofgem geocodes
#            installations to 2001 LSOAs, the Census and the IoD use 2011
#            LSOAs: the units are the connected components of the ONS
#            2001-2011 lookup (1:1 for 97% of LSOAs; split and merged LSOAs
#            are re-aggregated to the coarser of the two geographies)
#  Outcome:  number of domestic PV installations of up to 10 kW accredited
#            under the Feed-in Tariff (FiT) scheme, April 2010 - March 2019
#            (Ofgem Central FiT Register, installation report 30 June 2026).
#            Larger "domestic" installations (10-50 kW, 2.8% of those after
#            the 2016 reform) are landlord systems on blocks of flats, not
#            household adoptions; counts of all capacities and of up to 4 kW
#            are kept for the robustness checks
#  Covariates: 2011 Census (baseline, pre-scheme), English Indices of
#            Deprivation 2019, VOA Council Tax stock of properties 2015
#            (band and build period), DESNZ domestic electricity
#            consumption 2010 (before the scheme), roof form and floor area
#            from Energy Performance Certificates lodged before the reform,
#            conservation areas and listed buildings in force at the reform
#            (planning.data.gov.uk) and the leasehold share of houses sold
#            in 2008-2015 (HM Land Registry Price Paid)
#  Weights:  queen contiguity of the generalised (20 m) LSOA boundaries
#
#  Requires the duckdb command-line client to read the Ofgem workbooks
#  (1M rows; readxl takes several minutes on them) and the EPC and Land
#  Registry files, and an EPC API token (environment variable EPC_TOKEN or
#  the file ~/.epc_token) for the bulk EPC download (about 8 GB).
#  Writes data/solar/solar_lsoa.rds, W_lsoa.rds and lsoa_map.rds
# =============================================================================

library(sf)
library(spdep)
library(Matrix)
library(dplyr)
library(readr)
library(tidyr)

raw_dir <- "../data/solar/raw"
out_dir <- "../data/solar"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

download_once <- function(url, dest) {
  if (!file.exists(dest)) download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

# -----------------------------------------------------------------------------
# 1. Feed-in Tariff register: domestic PV installations per 2001 LSOA
# -----------------------------------------------------------------------------
ofgem <- "https://www.ofgem.gov.uk/sites/default/files/2026-07/feed-in_tariff_installation_report_"
fit_parts <- c(part_1 = "part_1_062026-20260716145808",
               part_2 = "part_2_062026-20260716145748_0",
               part_3 = "part_3_062026-20260716145832")
fit_xlsx <- file.path(raw_dir, paste0("fit_", names(fit_parts), ".xlsx"))
for (k in seq_along(fit_parts))
  download_once(paste0(ofgem, fit_parts[k], ".xlsx"), fit_xlsx[k])

fit_lsoa_csv <- file.path(raw_dir, "fit_pv_lsoa01_kw.csv")
fit_cov_csv  <- file.path(raw_dir, "fit_pv_coverage_kw.csv")
if (!file.exists(fit_lsoa_csv)) {
  # The header is on row 5 of each workbook; one row per accredited installation.
  # pv_dom and pv_dom_late count installations of up to 10 kW (the FiT band
  # boundary): the 10-50 kW "domestic" installations are mostly landlord
  # systems on blocks of flats, commissioned in batches
  sql <- sprintf("
    LOAD excel;
    CREATE TABLE fit AS
      SELECT * FROM read_xlsx('%s', header=true, range='A5:Z2000000', all_varchar=true)
      UNION ALL BY NAME SELECT * FROM read_xlsx('%s', header=true, range='A5:Z2000000', all_varchar=true)
      UNION ALL BY NAME SELECT * FROM read_xlsx('%s', header=true, range='A5:Z2000000', all_varchar=true);
    CREATE TABLE pv AS
      SELECT *, \"Installation Type\" = 'Domestic' AS dom,
             try_cast(\"Installed capacity\" AS DOUBLE) AS kw,
             try_strptime(\"Commissioning date\", '%%d/%%m/%%Y') >= DATE '2016-02-08' AS late
      FROM fit WHERE Technology = 'Photovoltaic' AND \"Extension (Y/N)\" = 'N';
    COPY (
      SELECT \"LLSOA Code\" AS lsoa01,
             count(*) FILTER (WHERE dom AND kw <= 10) AS pv_dom,
             sum(kw) FILTER (WHERE dom AND kw <= 10) AS pv_dom_kw,
             count(*) FILTER (WHERE dom AND kw <= 10 AND late) AS pv_dom_late,
             count(*) FILTER (WHERE dom) AS pv_dom_all,
             count(*) FILTER (WHERE dom AND late) AS pv_dom_late_all,
             count(*) FILTER (WHERE dom AND kw <= 4 AND late) AS pv_dom_late_4,
             count(*) FILTER (WHERE NOT dom) AS pv_nondom
      FROM pv WHERE \"LLSOA Code\" LIKE 'E01%%'
      GROUP BY 1
    ) TO '%s' (HEADER);
    COPY (
      SELECT count(*) AS n_dom_england,
             count(\"LLSOA Code\") AS n_dom_with_lsoa,
             count(*) FILTER (WHERE late) AS n_late_england,
             count(\"LLSOA Code\") FILTER (WHERE late) AS n_late_with_lsoa,
             count(*) FILTER (WHERE kw > 10) AS n_dom_gt10,
             count(\"LLSOA Code\") FILTER (WHERE late AND kw > 10) AS n_late_gt10_with_lsoa,
             count(*) FILTER (WHERE kw IS NULL) AS n_dom_no_kw
      FROM pv WHERE dom AND \"Installation Country\" = 'England'
    ) TO '%s' (HEADER);",
    fit_xlsx[1], fit_xlsx[2], fit_xlsx[3], fit_lsoa_csv, fit_cov_csv)
  stopifnot(system2("duckdb", c("-c", shQuote(sql))) == 0)
}
fit <- read_csv(fit_lsoa_csv, show_col_types = FALSE)
fit_vars <- c("pv_dom", "pv_dom_kw", "pv_dom_late", "pv_dom_all", "pv_dom_late_all",
              "pv_dom_late_4", "pv_nondom")

# -----------------------------------------------------------------------------
# 2. 2011 Census tables (Nomis API), all 2011 LSOAs (TYPE298), Total rows only
# -----------------------------------------------------------------------------
nomis_get <- function(id, cells = NULL, extra = "") {
  dest <- file.path(raw_dir, paste0("census2011_", id, ".csv"))
  if (!file.exists(dest)) {
    base <- sprintf(paste0("https://www.nomisweb.co.uk/api/v01/dataset/%s.data.csv?",
                           "geography=TYPE298&rural_urban=0&measures=20100",
                           "&select=geography_code,cell,cell_name,obs_value%s"), id, extra)
    if (!is.null(cells)) base <- paste0(base, "&cell=", paste(cells, collapse = ","))
    pages <- list(); off <- 0
    repeat {
      for (attempt in 1:5) {                 # the API drops connections now and then
        pg <- tryCatch(read_csv(paste0(base, "&recordoffset=", off), show_col_types = FALSE),
                       error = function(e) NULL)
        if (!is.null(pg)) break
        Sys.sleep(10)
      }
      stopifnot(!is.null(pg))
      if (nrow(pg) == 0) break
      pages[[length(pages) + 1]] <- pg
      off <- off + nrow(pg)
      if (nrow(pg) < 25000) break
    }
    write_csv(bind_rows(pages), dest)
  }
  read_csv(dest, show_col_types = FALSE) |>
    filter(startsWith(GEOGRAPHY_CODE, "E01"))
}
wide <- function(d, labels) {
  d |>
    filter(CELL_NAME %in% names(labels)) |>
    mutate(var = labels[CELL_NAME]) |>
    group_by(lsoa11 = GEOGRAPHY_CODE, var) |>
    summarise(v = sum(OBS_VALUE), .groups = "drop") |>
    pivot_wider(names_from = var, values_from = v)
}

ks401 <- wide(nomis_get("NM_618_1"), c(
  "Household spaces with at least one usual resident" = "households",
  "Household spaces with no usual residents" = "empty_spaces",
  "Whole house or bungalow: Detached" = "detached",
  "Whole house or bungalow: Semi-detached" = "semi",
  "Whole house or bungalow: Terraced (including end-terrace)" = "terraced",
  "Flat, maisonette or apartment: Purpose-built block of flats or tenement" = "flat",
  "Flat, maisonette or apartment: Part of a converted or shared house (including bed-sits)" = "flat",
  "Flat, maisonette or apartment: In a commercial building" = "flat",
  "Caravan or other mobile or temporary structure" = "caravan"))
ks402 <- wide(nomis_get("NM_619_1"), c(
  "All households" = "hh_tenure",
  "Owned" = "owned",
  "Shared ownership (part owned and part rented)" = "shared_own",
  "Social rented" = "social_rent",
  "Private rented" = "private_rent"))
ks501 <- wide(nomis_get("NM_623_1"), c(
  "All categories: Highest level of qualification" = "pop16",
  "Highest level of qualification: Level 4 qualifications and above" = "level4"))
ks102 <- wide(nomis_get("NM_145_1"), c(
  "All usual residents" = "pop",
  "Age 65 to 74" = "age65", "Age 75 to 84" = "age65",
  "Age 85 to 89" = "age65", "Age 90 and over" = "age65"))
qs415 <- wide(nomis_get("NM_547_1"), c(
  "All categories: Type of central heating in household" = "hh_heat",
  "Gas central heating" = "gas_heat"))
qs102 <- wide(nomis_get("NM_143_1"), c(
  "Area Hectares" = "area_ha"))
ks611 <- wide(nomis_get("NM_628_1", cells = c(0, 1, 4), extra = "&c_sex=0"), c(
  "All categories: NS-SeC" = "pop_nssec",
  "1. Higher managerial, administrative and professional occupations" = "nssec12",
  "2. Lower managerial, administrative and professional occupations" = "nssec12"))
ks105 <- wide(nomis_get("NM_605_1", cells = c(0, 5, 8, 10, 12)), c(
  "All categories: Household composition" = "hh_comp",
  "One family only: Married or same-sex civil partnership couple: Dependent children" = "hh_children",
  "One family only: Cohabiting couple: Dependent children" = "hh_children",
  "One family only: Lone parent: Dependent children" = "hh_children",
  "Other household types: With dependent children" = "hh_children"))
ks101 <- wide(nomis_get("NM_144_1", cells = 4), c(
  "Lives in a communal establishment" = "communal"))

# -----------------------------------------------------------------------------
# 2b. VOA Council Tax stock of properties, 2015 (2011 LSOAs): Council Tax band
#     (1991 values) and build period. Counts are rounded to 10; "-" marks
#     suppressed small counts, read as 0.
# -----------------------------------------------------------------------------
dir.create(file.path(raw_dir, "voa"), showWarnings = FALSE)
voa_num <- function(x) { x <- gsub(",", "", x); x[x %in% c("-", "N/A", "..")] <- "0"; as.numeric(x) }
voa_band <- read_csv(download_once(
  "https://assets.publishing.service.gov.uk/media/5a801553ed915d74e622c5d6/Table_CTSOP1.1_2015.csv",
  file.path(raw_dir, "voa", "ctsop1_1_2015.csv")), col_types = cols(.default = "c")) |>
  filter(GEOGRAPHY == "LSOA", startsWith(ECODE, "E01")) |>
  transmute(lsoa11 = ECODE,
            band_eh = voa_num(BAND_E) + voa_num(BAND_F) + voa_num(BAND_G) + voa_num(BAND_H),
            voa_props = voa_num(ALL_PROPERTIES))
voa_age <- read_csv(download_once(
  "https://assets.publishing.service.gov.uk/media/5a81b5c740f0b62305b906aa/Table_CTSOP4.1_2015.csv",
  file.path(raw_dir, "voa", "ctsop4_1_2015.csv")), col_types = cols(.default = "c")) |>
  filter(GEOGRAPHY == "LSOA", BAND == "All", startsWith(ECODE, "E01")) |>
  transmute(lsoa11 = ECODE, pre1919 = voa_num(BP_PRE_1900) + voa_num(BP_1900_1918),
            voa_props_bp = voa_num(ALL_PROPERTIES))

# -----------------------------------------------------------------------------
# 2c. DESNZ domestic electricity consumption by LSOA, 2010 (before the scheme;
#     later years net out self-consumption of the panels), on 2011 LSOAs
# -----------------------------------------------------------------------------
dir.create(file.path(raw_dir, "desnz"), showWarnings = FALSE)
elec_xlsx <- download_once(
  "https://assets.publishing.service.gov.uk/media/69427b7bd8156a816c419351/LSOA_domestic_elec_2010-2024.xlsx",
  file.path(raw_dir, "desnz", "lsoa_dom_elec_2010_2024.xlsx"))
elec_csv <- file.path(raw_dir, "desnz", "lsoa_dom_elec_2010.csv")
if (!file.exists(elec_csv)) {
  sql <- sprintf(paste("LOAD excel; COPY (SELECT * FROM read_xlsx('%s', sheet='2010',",
                       "header=true, range='A5:J60000', all_varchar=true)) TO '%s' (HEADER);"),
                 elec_xlsx, elec_csv)
  stopifnot(system2("duckdb", c("-c", shQuote(sql))) == 0)
}
elec11 <- read_csv(elec_csv, col_types = cols(.default = "c")) |>
  select(lsoa11 = 5, meters = 7, kwh = 8) |>
  filter(!is.na(lsoa11), startsWith(lsoa11, "E01")) |>
  transmute(lsoa11, elec_kwh = as.numeric(kwh) / as.numeric(meters))

# -----------------------------------------------------------------------------
# 3. English Indices of Deprivation 2019 (2011 LSOAs)
# -----------------------------------------------------------------------------
iod_csv <- download_once(
  "https://assets.publishing.service.gov.uk/media/5dc407b440f0b6379a7acc8d/File_7_-_All_IoD2019_Scores__Ranks__Deciles_and_Population_Denominators_3.csv",
  file.path(raw_dir, "iod", "iod2019_file7.csv"))
iod <- read_csv(iod_csv, show_col_types = FALSE) |>
  transmute(lsoa11 = `LSOA code (2011)`,
            lad19cd = `Local Authority District code (2019)`,
            lad19nm = `Local Authority District name (2019)`,
            income_dep = `Income Score (rate)`,
            iod_pop =`Total population: mid 2015 (excluding prisoners)`)

# -----------------------------------------------------------------------------
# 4. Geography: 2011 boundaries and population-weighted centroids, and the
#    2001-2011 LSOA lookup
# -----------------------------------------------------------------------------
arcgis <- "https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services"
arcgis_fetch <- function(svc, fields, code_field, dest) {
  if (!file.exists(dest)) {
    where <- utils::URLencode(sprintf("%s LIKE 'E01%%'", code_field), reserved = TRUE)
    n <- jsonlite::fromJSON(sprintf("%s/%s/FeatureServer/0/query?where=%s&returnCountOnly=true&f=json",
                                    arcgis, svc, where))$count
    pages <- lapply(seq(0, n - 1, by = 2000), function(off)
      st_read(sprintf(paste0("%s/%s/FeatureServer/0/query?where=%s&outFields=%s",
                             "&resultOffset=%d&resultRecordCount=2000&orderByFields=%s",
                             "&outSR=27700&f=geojson"),
                      arcgis, svc, where, fields, off, code_field), quiet = TRUE))
    st_write(do.call(rbind, pages), dest, delete_dsn = TRUE, quiet = TRUE)
  }
  st_read(dest, quiet = TRUE)
}
dir.create(file.path(raw_dir, "geo"), showWarnings = FALSE)
bgc <- arcgis_fetch("LSOA_Dec_2011_Boundaries_Generalised_Clipped_BGC_EW_V3",
                    "LSOA11CD,LSOA11NM", "LSOA11CD",
                    file.path(raw_dir, "geo", "lsoa11_bgc.gpkg"))
pwc <- arcgis_fetch("LSOA_Dec_2011_PWC_in_England_and_Wales_2022",
                    "lsoa11cd", "lsoa11cd",
                    file.path(raw_dir, "geo", "lsoa11_pwc.gpkg"))
st_crs(bgc) <- 27700; st_crs(pwc) <- 27700   # GeoJSON drops the BNG tag

pwc_xy <- st_coordinates(pwc)
pwc_df <- tibble(lsoa11 = pwc$lsoa11cd, east = pwc_xy[, 1], north = pwc_xy[, 2])

lk_csv <- file.path(raw_dir, "geo", "lsoa01_lsoa11_lookup.csv")
if (!file.exists(lk_csv)) {
  svc <- "LSOA01_LSOA11_LAD11_EW_LU_ddfe1cd1c2784c9b991cded95bc915a9"
  pages <- lapply(seq(0, 36000, by = 2000), function(off) jsonlite::fromJSON(sprintf(paste0(
    "%s/%s/FeatureServer/0/query?where=1%%3D1&outFields=LSOA01CD,LSOA11CD,CHGIND",
    "&resultOffset=%d&resultRecordCount=2000&orderByFields=FID&f=json"),
    arcgis, svc, off))$features$attributes)
  write_csv(bind_rows(pages), lk_csv)
}
lk <- read_csv(lk_csv, show_col_types = FALSE) |>
  filter(startsWith(LSOA11CD, "E01")) |>
  distinct(LSOA01CD, LSOA11CD, CHGIND)

# Four 2011 LSOAs created by 2001-2011 splits are absent from the ONS lookup:
# link each to the 2001 parents of the adjacent changed 2011 LSOAs of the same
# MSOA (this can only merge units, never drop installations)
miss <- setdiff(iod$lsoa11, lk$LSOA11CD)
if (length(miss)) {
  msoa <- function(cd) sub(".$", "", bgc$LSOA11NM[match(cd, bgc$LSOA11CD)])
  adj <- st_intersects(bgc[match(miss, bgc$LSOA11CD), ], bgc)
  add <- bind_rows(lapply(seq_along(miss), function(k) {
    nbr <- setdiff(bgc$LSOA11CD[adj[[k]]], miss[k])
    nbr <- nbr[msoa(nbr) == msoa(miss[k]) & nbr %in% lk$LSOA11CD[lk$CHGIND != "U"]]
    stopifnot(length(nbr) > 0)
    tibble(LSOA01CD = unique(lk$LSOA01CD[lk$LSOA11CD %in% nbr]),
           LSOA11CD = miss[k], CHGIND = "X")
  }))
  lk <- bind_rows(lk, add)
}

# -----------------------------------------------------------------------------
# 4b. Barriers to installation measured before the reform (8 Feb 2016), on
#     2011 LSOAs through the 2011-based National Statistics Postcode Lookup:
#     roof form and floor area from Energy Performance Certificates (latest
#     certificate per dwelling lodged before the reform), conservation areas
#     and listed buildings in force at the reform, and the leasehold share of
#     houses (latest sale per property in Price Paid, 2008-2015)
# -----------------------------------------------------------------------------
cut_date <- as.Date("2016-02-08")
duck <- function(sql) stopifnot(system2("duckdb", c("-c", shQuote(sql))) == 0)
nspl_csv <- file.path(raw_dir, "geo", "nspl11", "Data", "NSPL_AUG_2024_UK.csv")
if (!file.exists(nspl_csv)) {
  dir.create(dirname(dirname(nspl_csv)), showWarnings = FALSE)
  zip <- download_once("https://www.arcgis.com/sharing/rest/content/items/c5afedb9204a47e99559a4880feddcb1/data",
                       file.path(dirname(dirname(nspl_csv)), "nspl11_aug2024.zip"))
  stopifnot(system2("unzip", c("-q", "-o", shQuote(zip), "Data/NSPL_AUG_2024_UK.csv",
                               "-d", shQuote(dirname(dirname(nspl_csv))))) == 0)
}
nspl_sql <- sprintf("CREATE TABLE nspl AS SELECT upper(replace(pcds, ' ', '')) AS pc, lsoa11
  FROM read_csv('%s', all_varchar = true) WHERE lsoa11 LIKE 'E01%%';", nspl_csv)

# EPC: the bulk domestic file needs an API token of the EPC data service
epc_dir <- file.path(raw_dir, "epc"); dir.create(epc_dir, showWarnings = FALSE)
epc_zip <- file.path(epc_dir, "domestic-csv.zip")
if (!file.exists(epc_zip)) {
  tok <- Sys.getenv("EPC_TOKEN")
  if (!nzchar(tok) && file.exists("~/.epc_token"))
    tok <- gsub("[[:space:]]", "", readLines("~/.epc_token", warn = FALSE)[1])
  if (!nzchar(tok)) stop("the EPC bulk download needs an API token in EPC_TOKEN or ~/.epc_token")
  op <- options(timeout = 36000)
  download.file("https://api.get-energy-performance-data.communities.gov.uk/api/files/domestic/csv",
                epc_zip, mode = "wb", quiet = TRUE, headers = c(Authorization = paste("Bearer", tok)))
  options(op); rm(tok)
}
epc_files <- file.path(epc_dir, sprintf("certificates-%d.csv", 2008:2016))
for (f in epc_files[!file.exists(epc_files)])
  stopifnot(system2("unzip", c("-q", "-n", shQuote(epc_zip), basename(f), "-d", shQuote(epc_dir))) == 0)
epc_csv <- file.path(epc_dir, "epc_lsoa11_pre2016.csv")
epc_pc_csv <- file.path(epc_dir, "epc_postcode_dwellings_pre2016.csv")
if (!file.exists(epc_csv) || !file.exists(epc_pc_csv)) duck(sprintf("
  SET memory_limit = '16GB';
  CREATE TABLE epc AS
  SELECT uprn, upper(replace(postcode, ' ', '')) AS pc, lower(address) AS addr, lodgement_date,
         roof_description AS rd, try_cast(total_floor_area AS DOUBLE) AS tfa
  FROM read_csv([%s], all_varchar = true, union_by_name = true)
  WHERE country = 'England' AND lodgement_date < '%s';
  CREATE TABLE last AS SELECT * EXCLUDE (rn) FROM (
    SELECT *, row_number() OVER (PARTITION BY coalesce(uprn, pc || '|' || addr)
                                 ORDER BY lodgement_date DESC) AS rn FROM epc) WHERE rn = 1;
  %s
  COPY (
    SELECT n.lsoa11, count(*) AS n_epc,
      sum(CASE WHEN rd ILIKE '%%dwelling above%%' OR rd ILIKE '%%premises above%%' THEN 1 ELSE 0 END) AS no_roof,
      sum(CASE WHEN rd ILIKE 'Flat%%' THEN 1 ELSE 0 END) AS flat_roof,
      sum(CASE WHEN rd ILIKE 'Roof room%%' THEN 1 ELSE 0 END) AS roof_room,
      sum(CASE WHEN tfa > 0 THEN tfa END) AS tfa_sum, sum(CASE WHEN tfa > 0 THEN 1 ELSE 0 END) AS tfa_n
    FROM last l JOIN nspl n USING (pc) GROUP BY n.lsoa11
  ) TO '%s' (HEADER);
  COPY (SELECT pc, count(*) AS n_dwell FROM last GROUP BY pc) TO '%s' (HEADER);",
  paste0("'", epc_files, "'", collapse = ", "), cut_date, nspl_sql, epc_csv, epc_pc_csv))
epc11 <- read_csv(epc_csv, show_col_types = FALSE) |>
  select(lsoa11, n_epc, no_roof, flat_roof, roof_room, tfa_sum, tfa_n)

# Conservation areas designated before the reform and not ended by then: share
# of the EPC dwellings, located at their postcode centroid, that lie inside one
plan_dir <- file.path(raw_dir, "planning"); dir.create(plan_dir, showWarnings = FALSE)
ca <- st_read(download_once("https://files.planning.data.gov.uk/dataset/conservation-area.geojson",
                            file.path(plan_dir, "conservation-area.geojson")), quiet = TRUE)
in_force <- function(start, end) {
  start <- as.Date(ifelse(nzchar(start), start, NA)); end <- as.Date(ifelse(nzchar(end), end, NA))
  (is.na(start) | start < cut_date) & (is.na(end) | end >= cut_date)
}
ca <- st_make_valid(st_transform(ca[in_force(ca$designation.date, ca$end.date), ], 27700))
nspl_xy <- read_csv(nspl_csv, col_select = c(pcds, oseast1m, osnrth1m, lsoa11), show_col_types = FALSE) |>
  filter(startsWith(lsoa11, "E01"), !is.na(oseast1m)) |>
  transmute(pc = toupper(gsub(" ", "", pcds)), oseast1m, osnrth1m, lsoa11)
pts <- read_csv(epc_pc_csv, show_col_types = FALSE) |> inner_join(nspl_xy, by = "pc")
pts$in_ca <- lengths(st_intersects(st_as_sf(pts, coords = c("oseast1m", "osnrth1m"), crs = 27700), ca)) > 0
ca11 <- pts |> group_by(lsoa11) |> summarise(dw = sum(n_dwell), dw_ca = sum(n_dwell * in_ca), .groups = "drop")

# Listed buildings listed before the reform and not removed by then
lb <- read_csv(download_once("https://files.planning.data.gov.uk/dataset/listed-building.csv",
                             file.path(plan_dir, "listed-building.csv")),
               col_select = c(point, `start-date`, `end-date`), col_types = cols(.default = "c")) |>
  filter(!is.na(point), nzchar(point), in_force(coalesce(`start-date`, ""), coalesce(`end-date`, "")))
hit <- st_intersects(st_transform(st_as_sf(lb, wkt = "point", crs = 4326), 27700), bgc)
lb11 <- tibble(lsoa11 = bgc$LSOA11CD[unlist(hit[lengths(hit) == 1])]) |> count(lsoa11, name = "listed")

# Leasehold: latest standard (category A) sale of each property before the
# reform, 2008-2015; share of houses (detached, semi-detached, terraced) held
# on a lease
lr_dir <- file.path(raw_dir, "landreg"); dir.create(lr_dir, showWarnings = FALSE)
for (yr in 2008:2015)
  download_once(sprintf("https://price-paid-data.publicdata.landregistry.gov.uk/pp-%d.csv", yr),
                file.path(lr_dir, sprintf("pp-%d.csv", yr)))
lease_csv <- file.path(lr_dir, "lease_lsoa11_pre2016.csv")
if (!file.exists(lease_csv)) duck(sprintf("
  CREATE TABLE pp AS SELECT * FROM read_csv('%s/pp-20*.csv', header = false, all_varchar = true,
    names = ['tid','price','date','postcode','ptype','new','duration','paon','saon','street','locality','town','district','county','cat','status']);
  CREATE TABLE prop AS SELECT * EXCLUDE (rn) FROM (
    SELECT upper(replace(postcode, ' ', '')) AS pc, ptype, duration,
           row_number() OVER (PARTITION BY postcode, paon, saon, street ORDER BY date DESC) AS rn
    FROM pp WHERE cat = 'A' AND postcode IS NOT NULL AND date < '%s') WHERE rn = 1;
  %s
  COPY (SELECT lsoa11, sum(CASE WHEN ptype IN ('D','S','T') THEN 1 ELSE 0 END) AS pp_houses,
               sum(CASE WHEN ptype IN ('D','S','T') AND duration = 'L' THEN 1 ELSE 0 END) AS pp_lease_houses
        FROM prop JOIN nspl USING (pc) GROUP BY lsoa11) TO '%s' (HEADER);",
  lr_dir, cut_date, nspl_sql, lease_csv))
lease11 <- read_csv(lease_csv, show_col_types = FALSE)
if (!"pp_houses" %in% names(lease11))      # file written by an earlier version
  lease11 <- rename(lease11, pp_houses = houses, pp_lease_houses = lease_houses)
lease11 <- select(lease11, lsoa11, pp_houses, pp_lease_houses)

barriers11 <- iod |> select(lsoa11) |>
  left_join(epc11, by = "lsoa11") |> left_join(ca11, by = "lsoa11") |>
  left_join(lb11, by = "lsoa11") |> left_join(lease11, by = "lsoa11") |>
  mutate(across(-lsoa11, ~ coalesce(.x, 0)))
barrier_vars <- setdiff(names(barriers11), "lsoa11")

# -----------------------------------------------------------------------------
# 5. Harmonised units: connected components of the 2001-2011 lookup graph
# -----------------------------------------------------------------------------
nodes <- c(paste0("a", unique(lk$LSOA01CD)), paste0("b", unique(lk$LSOA11CD)))
i01 <- match(paste0("a", lk$LSOA01CD), nodes)
i11 <- match(paste0("b", lk$LSOA11CD), nodes)
lab <- seq_along(nodes)
repeat {                       # label propagation: each node takes the minimum
  m <- pmin(lab[i01], lab[i11])
  new <- lab
  for (ix in list(i01, i11)) {
    lm <- tapply(m, ix, min)
    k  <- as.integer(names(lm))
    new[k] <- pmin(new[k], lm)
  }
  if (identical(new, lab)) break
  lab <- new
}
lk$comp <- lab[i01]
units <- lk |>
  group_by(comp) |>
  summarise(n01 = n_distinct(LSOA01CD), n11 = n_distinct(LSOA11CD),
            unit = if (n11 == 1) first(LSOA11CD) else if (n01 == 1) first(LSOA01CD)
                   else paste0("X", min(LSOA11CD)),
            .groups = "drop")
lk <- left_join(lk, select(units, comp, unit), by = "comp")
unit_of01 <- distinct(lk, LSOA01CD, unit)
unit_of11 <- distinct(lk, LSOA11CD, unit)

fit <- fit |> left_join(unit_of01, by = c("lsoa01" = "LSOA01CD"))
fit_unmatched <- filter(fit, is.na(unit))
fit_u <- fit |>
  filter(!is.na(unit)) |>
  group_by(unit) |>
  summarise(across(all_of(fit_vars), sum), .groups = "drop")

# -----------------------------------------------------------------------------
# 6. Analysis dataset: Census counts summed over the 2011 LSOAs of each unit
# -----------------------------------------------------------------------------
l11 <- iod |>
  left_join(ks401, by = "lsoa11") |>
  left_join(ks402, by = "lsoa11") |>
  left_join(ks501, by = "lsoa11") |>
  left_join(ks102, by = "lsoa11") |>
  left_join(qs415, by = "lsoa11") |>
  left_join(qs102, by = "lsoa11") |>
  left_join(ks611, by = "lsoa11") |>
  left_join(ks105, by = "lsoa11") |>
  left_join(ks101, by = "lsoa11") |>
  left_join(voa_band, by = "lsoa11") |>
  left_join(voa_age, by = "lsoa11") |>
  left_join(elec11, by = "lsoa11") |>
  left_join(pwc_df, by = "lsoa11") |>
  left_join(barriers11, by = "lsoa11") |>
  left_join(unit_of11, by = c("lsoa11" = "LSOA11CD"))
stopifnot(nrow(l11) == 32844, !anyNA(l11$unit), !anyNA(l11$households),
          !anyNA(select(l11, nssec12, hh_children, communal, band_eh, pre1919, elec_kwh)))

solar <- l11 |>
  group_by(unit) |>
  summarise(n_lsoa11 = n(), lad19cd = first(lad19cd), lad19nm = first(lad19nm),
            east = weighted.mean(east, pop), north = weighted.mean(north, pop),
            income_dep = sum(income_dep * iod_pop) / sum(iod_pop),
            elec_kwh = weighted.mean(elec_kwh, households),
            across(c(households, detached, semi, terraced, flat, caravan,
                     hh_tenure, owned, shared_own, social_rent, private_rent,
                     pop16, level4, pop, age65, hh_heat, gas_heat, area_ha,
                     empty_spaces, pop_nssec, nssec12, hh_comp, hh_children,
                     communal, band_eh, voa_props, pre1919, voa_props_bp,
                     all_of(barrier_vars)), sum),
            .groups = "drop") |>
  left_join(select(units, unit, n01), by = "unit") |>
  left_join(fit_u, by = "unit") |>
  mutate(across(all_of(fit_vars), ~ coalesce(.x, 0)),
         houses       = detached + semi + terraced,
         spaces       = houses + flat + caravan,   # household spaces, incl. vacant
         sh_owner     = (owned + shared_own) / hh_tenure,
         sh_social    = social_rent / hh_tenure,
         sh_private   = private_rent / hh_tenure,
         sh_detached  = detached / spaces,
         sh_flat      = flat / spaces,
         sh_level4    = level4 / pop16,
         sh_age65     = age65 / pop,
         sh_offgas    = 1 - gas_heat / hh_heat,
         density      = pop / (area_ha / 100),   # residents per km^2
         sh_nssec12   = nssec12 / pop_nssec,
         sh_children  = hh_children / hh_comp,
         sh_communal  = communal / pop,
         sh_empty     = empty_spaces / (households + empty_spaces),
         sh_band_eh   = band_eh / voa_props,
         sh_pre1919   = pre1919 / voa_props_bp,
         # barriers before the reform; NA where a unit has no EPC dwelling or
         # no house sale to measure them on
         epc_noroof   = ifelse(n_epc > 0, no_roof / n_epc, NA),      # dwelling or premises above
         epc_flatroof = ifelse(n_epc > 0, flat_roof / n_epc, NA),
         epc_roofroom = ifelse(n_epc > 0, roof_room / n_epc, NA),
         epc_floor    = ifelse(tfa_n > 0, log(tfa_sum / tfa_n), NA), # log mean floor area (m^2)
         epc_cover    = ifelse(n_epc > 0, n_epc / households, NA),   # EPC dwellings per household
         ca_share     = ifelse(dw > 0, dw_ca / dw, NA),              # EPC dwellings in a conservation area
         listed_rate  = 1000 * listed / households,                  # listed buildings per 1,000 households
         lease_house  = ifelse(pp_houses > 0, pp_lease_houses / pp_houses, NA),
         pv_rate      = 100 * pv_dom / households) |>
  arrange(unit)

# -----------------------------------------------------------------------------
# 7. Spatial weights: queen contiguity, islands linked to nearest centroid
# -----------------------------------------------------------------------------
bgc$unit <- unit_of11$unit[match(bgc$LSOA11CD, unit_of11$LSOA11CD)]
multi <- bgc$unit %in% solar$unit[solar$n_lsoa11 > 1]
geo <- rbind(
  select(bgc[!multi, ], unit),
  bgc[multi, ] |> group_by(unit) |> summarise(.groups = "drop"))
geo <- geo[match(solar$unit, geo$unit), ]
stopifnot(identical(geo$unit, solar$unit))

nb <- poly2nb(geo, queen = TRUE, snap = 1)
isl <- which(card(nb) == 0)
if (length(isl)) {
  knn1 <- knn2nb(knearneigh(cbind(solar$east, solar$north), k = 1))
  for (i in isl) nb[[i]] <- knn1[[i]]
  nb <- make.sym.nb(nb)
}
# row-standardised sparse W from a neighbour list
nb_to_W <- function(nb) {
  k <- card(nb)
  sparseMatrix(i = rep(seq_along(nb), k), j = unlist(nb), x = rep(1 / k, k),
               dims = c(length(nb), length(nb)))
}
W_queen <- nb_to_W(nb)

# k-nearest-neighbour alternatives for the sensitivity analysis
xy <- cbind(solar$east, solar$north)
W_knn <- lapply(c(k6 = 6, k10 = 10), function(k)
  nb_to_W(knn2nb(knearneigh(xy, k = k))))

saveRDS(solar, file.path(out_dir, "solar_lsoa.rds"))
saveRDS(list(queen = W_queen, nb_queen = nb, k6 = W_knn$k6, k10 = W_knn$k10,
             islands = solar$unit[isl]),
        file.path(out_dir, "W_lsoa.rds"))
saveRDS(st_simplify(geo, dTolerance = 50), file.path(out_dir, "lsoa_map.rds"))
saveRDS(unit_of11, file.path(out_dir, "unit_lookup.rds"))   # 2011 LSOA -> unit, for crime_data.r

cov <- read_csv(fit_cov_csv, show_col_types = FALSE)
cat("Domestic PV in England:", cov$n_dom_england, " with LSOA code:", cov$n_dom_with_lsoa,
    " assigned to units:", sum(solar$pv_dom_all), "(", sum(solar$pv_dom), "up to 10 kW )",
    " unmatched codes:", nrow(fit_unmatched), "(", sum(fit_unmatched$pv_dom_all), "installations )\n")
cat("After the reform:", cov$n_late_england, " with LSOA code:", cov$n_late_with_lsoa,
    sprintf("(%.1f%% without)", 100 * (1 - cov$n_late_with_lsoa / cov$n_late_england)),
    " over 10 kW among them:", cov$n_late_gt10_with_lsoa,
    " in units: all", sum(solar$pv_dom_late_all), " up to 10 kW", sum(solar$pv_dom_late),
    " up to 4 kW", sum(solar$pv_dom_late_4), "\n")
cat("Units without a pre-reform EPC:", sum(is.na(solar$epc_noroof)),
    " without a house sale:", sum(is.na(solar$lease_house)),
    " without a located EPC dwelling:", sum(is.na(solar$ca_share)), "\n")
cat("Units:", nrow(solar), " (2011 LSOAs:", sum(solar$n_lsoa11),
    "; units of >1 2011 LSOA:", sum(solar$n_lsoa11 > 1),
    "; of >1 2001 LSOA:", sum(solar$n01 > 1), ")\n")
cat("zeros:", sum(solar$pv_dom == 0), sprintf("(%.2f%%)", 100 * mean(solar$pv_dom == 0)),
    " ones:", sum(solar$pv_dom == 1),
    " mean:", round(mean(solar$pv_dom), 2), " var:", round(var(solar$pv_dom), 1), "\n")
cat("Queen neighbours: mean", round(mean(card(nb)), 2),
    " islands relinked:", length(isl), "\n")
