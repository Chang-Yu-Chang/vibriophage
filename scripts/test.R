
library(tidyverse)
source(here::here("scripts/simulate_models.R"))

# Shared baseline that keeps plankton alive (important!)
p_base <- list(
    # Human FOI (saturating)
    beta=1.2e-6, h_V=1e-5,
    # Vibrio & phage
    sigma=50, delta_V=0.30, tau=80, delta_P=0.20,
    mu = 3e-6,
    # Plankton (balanced so X,Z don’t crash)
    lambda=0.18, c_lambda=1.2, epsilon=0.50,
    eta=0.6, alpha=0.08, delta_Z=0.04, r=1.0, K=3.0, kappa_Z=10,
    # Infectivity mixing
    theta=0.8, g=0.0  # <— VZ immune in M3
)
i_M0 <- list(S=9.99e5, I=3, VF=60, VZ=0, P=0,   ZF=0,  ZV=0,   X=0)
i_M1 <- list(S=9.99e5, I=5, VF=80, VZ=3, P=0,   ZF=6,  ZV=0.3, X=1.5)
i_M2 <- list(S=9.99e5, I=5, VF=80, VZ=0, P=4e5, ZF=0,  ZV=0,   X=0)
i_M3 <- list(S=9.99e5, I=5, VF=80, VZ=3, P=4e5, ZF=6,  ZV=0.3, X=1.5)

# Choose the “strong phage, moderate detachment” you liked:
#mu_use <- 3e-6   # strong adsorption
eps_use <- 0.5   # ~2 d residence
times <- 0:120

runs <- tibble(model = c("M0","M1","M2","M3")) %>%
    mutate(init = list(i_M0, i_M1, i_M2, i_M3),
           params = list(p_base, p_base, p_base, p_base)) %>%
    # mutate(params = list(
    #     modifyList(p_base, list(mu = 0.0)),               # M0: no phage
    #     modifyList(p_base, list(mu = 0.0)),               # M1: no phage
    #     modifyList(p_base, list(mu = mu_use)),            # M2: phage-only
    #     modifyList(p_base, list(mu = mu_use, g = 0.0))    # M3: phage + immune VZ
    # ),
    mutate(sim = map2(model, 1:n(), ~{
        simulate_model(.x, times=times, params=params[[.y]], init=init[[.y]]) %>%
            mutate(V_tot = VF + ifelse("VZ" %in% names(.), VZ, 0)) %>%
            mutate(across(-model, as.numeric)) %>%
            dplyr::select(-model)
    })) %>%
    select(model, sim) %>%
    unnest(sim)

# Plots (peak too high in M0/M1; wiped in M2; “reasonable” in M3)
ggplot(runs, aes(time, I, color = model)) +
    geom_line() + facet_wrap(~model, scales="free_y") + theme_bw()
