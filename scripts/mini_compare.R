# mini_compare.R
# Synthetic test: generate daily cases from M4, then score M0–M4.

library(tidyverse)
library(lubridate)
source(here::here("scripts/simulate_models.R"))
source(here::here("scripts/observe_and_score.R"))

# Observation
# haiti_cholera <- read_csv(here::here("data/haiti_cholera.csv"))
# obs <- haiti_cholera %>%
#     filter(brief_desc == "Haiti_MSPP_NASA")

set.seed(1)
# Simulate 180 days
times <- 0:180
# Truth: M4 with a 3-stage protection
sim_truth <- simulate_model("M4", times = times, k = 3)
ex_truth  <- expected_daily_cases(sim_truth)

# Generate synthetic observed daily cases via NegBin (phi fixed)
phi <- 10
rho_true <- 0.2
obs_daily <- tibble(
    date  = as.Date("2010-10-17") + times,
    cases = rnbinom(length(times), mu = rho_true * ex_truth$mu, size = phi)
)

# Compare models
models <- c("M0","M1","M2","M3","M4")
scores <- compare_models_daily(models, simulate_model,
                               obs_daily_df = obs_daily,
                               start_date   = min(obs_daily$date),
                               times        = times,
                               k = 3, phi = phi)

print(scores %>% arrange(desc(NB_logscore)))
