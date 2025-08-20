# plot_simulations.R

library(tidyverse)
source(here::here("scripts/simulate_models.R"))

params_updated <- list(
    # Human infection scale (saturating FOI)
    beta   = 1.5e-6,     # conc^-1 day^-1  (keeps early hazard realistic)
    h_V    = 1e-5,       # half-saturation (1/h_V ~ 1e5 conc units)

    # Vibrio baseline
    sigma  = 50,         # shedding from I (keep as before)
    delta_V= 0.30,       # day^-1  (washout/mortality ~ 3.3 d)

    # Phage life history (moderate—not instant kill, not vanishing)
    mu     = 3e-7,       # L · virion^-1 · day^-1
    tau    = 80,         # burst size
    delta_P= 0.20,       # day^-1  (half-life ~3.5 d)

    # Plankton coupling (let VZ form & matter)
    lambda   = 1,    # VF–ZF encounter
    c_lambda = 10,     # colonization efficiency
    epsilon  = 0.5,    # VZ mean residence ~ 2 d
    theta    = 10,       # VZ contributes strongly to FOI
    g        = 1e-2,    # partial susceptibility of VZ to phage (M3)

    # Zooplankton–phyto (avoid X crash)
    alpha   = 0.08,     # LOWER grazing rate
    r       = 1.0,      # HIGHER phytoplankton growth
    K       = 3.0,      # BIGGER carrying capacity
    delta_Z = 0.04,     # LOWER zooplankton mortality
    eta     = 0.6,
    kappa_Z = 10,

    # Immunity loss left at your previous default
    omega  = 9.2e-4
)

set.seed(1)
times <- 0:100

tb <- tibble(model = paste0("M", 0:3)) %>%
    mutate(sim = purrr::map(model, ~simulate_model(.x, times = times, params = params_updated) %>% mutate(across(-model, as.numeric))
    )) %>%
    select(sim) %>% unnest(sim) %>%
    mutate(V = VF + VZ, Z = ZF + ZV)


p1 <- tb %>%
    ggplot(aes(time, ZF, color = model)) +
    geom_line() +
    facet_grid(~model, scales = "free_y") +
    theme_bw()

p2 <- tb %>%
    ggplot(aes(time, ZV, color = model)) +
    geom_line() +
    facet_grid(~model, scales = "free_y") +
    theme_bw()

p3 <- tb %>%
    ggplot(aes(time, VF, color = model)) +
    geom_line() +
    facet_grid(~model, scales = "free_y") +
    theme_bw()

p4 <- tb %>%
    ggplot(aes(time, VZ, color = model)) +
    geom_line() +
    facet_grid(~model, scales = "free_y") +
    theme_bw()

p5 <- tb %>%
    ggplot(aes(time, I, color = model)) +
    geom_line() +
    facet_grid(~model, scales = "free_y") +
    theme_bw()

cowplot::plot_grid(p1, p2, p3, p4, p5, ncol = 1, align = "v", axis = "rl")



# tb %>%
#     select(sim) %>% unnest(sim) %>%
#     ggplot(aes(time, I, color = model)) +
#     geom_line() +
#     facet_grid(~model, scales = "free_y") +
#     theme_bw()

#
# times <- 0:100
# m0 <- simulate_model("M0", times = times)
# m1 <- simulate_model("M1", times = times)
# m2 <- simulate_model("M2", times = times)
# m3 <- simulate_model("M3", times = times)
#
