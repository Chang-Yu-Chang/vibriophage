# obserbe_and_score.R
# Utilities to map states -> daily incidence, fit rho, and score models.

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----- Expected daily incidence from a simulation -----
# Uses the same saturating FOI as the ODEs:
#   lambda_H = beta * (VF + theta * VZ) / (1 + h_V * (VF + theta * VZ))
# Expected new infections per day = lambda_H * S
expected_daily_cases <- function(sim_df,
                                 beta  = attr(sim_df, "params")$beta %||% stop("beta missing"),
                                 theta = attr(sim_df, "params")$theta %||% 0,
                                 h_V   = attr(sim_df, "params")$h_V   %||% 0) {
    VZ <- if ("VZ" %in% names(sim_df)) sim_df$VZ else 0
    Veff <- sim_df$VF + theta * VZ
    lambda_H <- beta * Veff / (1 + h_V * Veff)
    mu <- pmax(lambda_H * sim_df$S, 0)  # expected daily cases (unreported)
    tibble(time = sim_df$time, mu = mu)
}

# ----- Fit a single reporting fraction rho (Poisson MLE) on aligned daily data -----
fit_rho_daily <- function(obs_daily_df, sim_df, start_date) {
    ex <- expected_daily_cases(sim_df) %>% mutate(date = as.Date(start_date) + time)
    cmp <- obs_daily_df %>% left_join(ex, by = "date") %>% drop_na()
    if (!nrow(cmp)) stop("No overlapping dates between observed and simulation.")
    sum(cmp$cases) / sum(cmp$mu)
}

# ----- Scores -----
rmse_daily <- function(y, mu_scaled) sqrt(mean((y - mu_scaled)^2, na.rm = TRUE))

# NB log-score using 'mu' parameterization; larger (less negative) is better
nb_logscore_daily <- function(y, mu_scaled, phi = 10) {
    mu_scaled <- pmax(mu_scaled, 1e-9)
    sum(dnbinom(y, mu = mu_scaled, size = phi, log = TRUE))
}

# ----- Convenience: compare a list of models against the same daily series -----
# simulate_model_fn must be your simulate_model(model, times, params, init)
compare_models_daily <- function(models, simulate_model_fn,
                                 obs_daily_df, start_date, times,
                                 phi = 10, params = list(), init = list()) {
    purrr::map_dfr(models, function(m) {
        sim <- simulate_model_fn(m, times = times, params = params, init = init)
        rho <- fit_rho_daily(obs_daily_df, sim, start_date)
        ex  <- expected_daily_cases(sim) %>% mutate(date = as.Date(start_date) + time)
        cmp <- obs_daily_df %>% left_join(ex, by = "date") %>% drop_na()
        tibble(
            model = m,
            rho_hat = rho,
            RMSE = rmse_daily(cmp$cases, rho * cmp$mu),
            NB_logscore = nb_logscore_daily(cmp$cases, rho * cmp$mu, phi = phi)
        )
    })
}
