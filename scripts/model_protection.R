# model_protection.R
# Two tiny models (No-Protection vs Protection) + B1 (time series) + B2 (robustness bars)

suppressPackageStartupMessages({
    library(deSolve)
    library(tidyverse)
    library(patchwork)
})

# -------------------- Output helper --------------------
expected_daily <- function(sim_df) {
    p  <- attr(sim_df, "params")
    df <- as_tibble(sim_df) %>% mutate(across(everything(), as.numeric))
    Veff     <- df$VF + df$VI + p[["theta"]]*df$VZ
    lambda_H <- p[["beta"]] * Veff
    tibble(time = df$time, daily = lambda_H * df$S)
}

# -------------------- Core ODEs --------------------
rhs_factory <- function() {
    function(t, y, p) {
        with(as.list(c(y, p)), {

            # Human compartments
            N <- S + I + R
            Veff <- VF + VI + theta*VZ # effective vibrio population size
            lambda_H <- beta * Veff # type I response

            dS <- b*N + omega*R - delta*S - lambda_H*S
            dI <- lambda_H*S - (gamma + delta_I + delta)*I
            dR <- gamma*I - (omega + delta)*R

            # Free Vibrio
            dVF <- sigma*I + kappa*VI - lambda*VF*Z0 - mu*VF*P - delta_V*VF

            # Zooplankton associated Vibrio
            dVZ <- lambda*VF*Z0 - epsilon*VZ*Z0 - delta_V*VZ

            # Immune Vibrio (rho = detachment fraction, g = immunity)
            dVI <- epsilon*VZ - kappa*VI - (1-g)*mu*VI*P - delta_V*VI

            # Phage
            dP <- tau*mu*(VF + (1-g)*VI)*P - delta_P*P

            list(c(dS, dI, dR, dVF, dVZ, dVI, dP))
        })
    }
}

# -------------------- Simulation --------------------
simulate_min <- function(times = 0:500, params = list(), init = list()) {
    p <- modifyList(list(
        # demographics (kept neutral for clarity)
        #b=6.85e-5, delta=3.8e-5, omega=1e-4,
        b=0, delta=0, omega=0,
        # infection
        theta=0,                # VZ does not directly infect here
        beta=1.2e-7,
        gamma=0.6,
        delta_I=0.05,
        # Vibrio & phage
        sigma=50, delta_V=0.2,
        mu=8e-5, tau=200, delta_P=0.02,  # strong phage pressure: key for contrast
        # Z coupling
        lambda=20,
        Z0=1,
        epsilon=1,            # moderate release: too small → very slow; too big → washed out
        # protection strength & window
        g=1,                    # immunity strength (1=full, 0=none)
        kappa=.1               # VI → VF return (~1/day; shorter if you want faster cycles)
    ), params)

    p <- unlist(p)   # <-- key step: force numerics

    y0 <- modifyList(list(
        S=1e6, I=0, R=0,
        VF=0, VZ=1e2, VI=0, P=4e5
    ), init)

    y <- unlist(y0)

    rhs <- rhs_factory()
    out <- ode(y=y, times=times, func=rhs, parms=p, method="lsoda") %>%
        as_tibble() %>%
        mutate(across(everything(), as.numeric))
    attr(out, "params") <- p
    out
}



# -------------------- Run scenarios --------------------
times = 0:500
mu = 1e-5
sim_no  <- simulate_min(times=times, params = list(mu=mu, g=0))
sim_yes <- simulate_min(times=times, params = list(mu=mu, g=1, kappa=log(2)/(5/24)))

# bind_rows(
#     sim_no  %>% mutate(model="No protection"),
#     sim_yes %>% mutate(model="Protection")
# ) %>%
#     pivot_longer(-c(model, time)) %>%
#     mutate(name = factor(name, c("S", "I", "R", "VF", "VI", "VZ", "P"))) %>%
#     ggplot() +
#     geom_line(aes(time, value, color=model)) +
#     facet_wrap(~name, scales = "free_y") +
#     coord_cartesian(clip = "off") +
#     theme_bw() +
#     theme() +
#     guides() +
#     labs()

# (A) Daily infections (main panel)
bind_rows(
    expected_daily(sim_no)  |> mutate(model="No protection (g=0)"),
    expected_daily(sim_yes) |> mutate(model="Protection (g=1)")
) |>
    ggplot(aes(time, daily, color=model)) +
    geom_line(linewidth=0.9, alpha = .8) +
    scale_y_continuous(labels=scales::comma) +
    labs(x="Time (days)", y="Expected daily infections",
         title="Immunity is necessary for recurring outbreaks in presence of phage",
         subtitle="g=0: phage purge free Vibrio → single small transient; g=1: VI refuge sustains cycles") +
    theme_bw() + theme(legend.position="top")
