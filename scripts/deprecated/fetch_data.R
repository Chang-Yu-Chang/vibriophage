# ============================================================
# Fetch WHO Global Cholera weekly data (ArcGIS FeatureServer)
# ============================================================

library(httr)
library(jsonlite)
library(dplyr)

# Example endpoint (replace with the actual FeatureServer layer URL)
# You’ll get this by inspecting the dashboard network calls in devtools.
# It usually looks like:
# "https://services.arcgis.com/.../FeatureServer/0/query"

url <- "https://example.com/arcgis/rest/services/cholera_weekly/FeatureServer/0/query"

# Query parameters: return all fields and all rows
params <- list(
    where = "1=1",              # no filter
    outFields = "*",            # get all columns
    f = "json",                 # format = json
    returnGeometry = "false"    # no map geometry, just attributes
)

# Fetch data
resp <- GET(url, query = params)
dat_json <- content(resp, as = "text", encoding = "UTF-8")

# Parse JSON into dataframe
dat <- fromJSON(dat_json, flatten = TRUE)$features
cholera_df <- dat$attributes %>% as_tibble()

# Inspect structure
glimpse(cholera_df)

# Example: rename time fields if available
cholera_df <- cholera_df %>%
    mutate(
        epi_week = as.integer(epi_week),   # some layers use epi_week or week_num
        report_date = as.Date(report_date/1000, origin="1970-01-01") # if timestamp in ms
    )

# Save locally if you want
# write.csv(cholera_df, "cholera_weekly.csv", row.names = FALSE)
