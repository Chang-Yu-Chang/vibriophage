# model_protection.R
# Two tiny models (No-Protection vs Protection) + B1 (time series) + B2 (robustness bars)

suppressPackageStartupMessages({
    library(deSolve)
    library(tidyverse)
    library(patchwork)
})

# -------------------- Output helper --------------------
expected_daily <- function(sim_df) {
    p <- attr(sim_df, "params")  # named numeric vector
    sim_tb <- as_tibble(sim_df) %>%
        mutate(across(everything(), as.numeric))  # strip attributes & classes
    Veff <- sim_tb$VF + p[["theta"]] * (sim_tb$VZ + sim_tb$VI)
    lambda_H <- p[["beta"]] * Veff / (1 + p[["h_V"]] * Veff)

    tibble(time = sim_tb$time, daily = lambda_H * sim_tb$S)
}
# -------------------- Core ODEs --------------------

rhs_factory <- function() {
    function(t, y, p) {
        with(as.list(c(y, p)), {

            # Human compartments
            N <- S + I + R
            Veff <- (VF + VI) + theta*VZ # effective vibrio population size
            lambda_H <- beta * Veff^2 / (1 + h_V * Veff^2) # type III infection

            dS <- b*N + omega*R - delta*S - lambda_H*S
            dI <- lambda_H*S - (gamma + delta_I + delta)*I
            dR <- gamma*I - (omega + delta)*R

            # Free Vibrio
            dVF <- sigma*I - delta_V*VF - lambda*VF*Z0 + (1-rho)*epsilon*VZ + kappa*VI - mu*VF*P

            # Zooplankton associated Vibrio
            dVZ <- lambda*VF*Z0 - epsilon*VZ - delta_V*VZ # - (1-g)*mu*VZ*P

            # Immune Vibrio (rho = detachment fraction, g = immunity)
            dVI <- rho*epsilon*VZ - kappa*VI - delta_V*VI - (1-g)*mu*VI*P

            # Phage
            dP <- tau*mu*(VF + (1-g)*VI + (1-g)*VZ)*P - delta_P*P

            list(c(dS, dI, dR, dVF, dVZ, dVI, dP))
        })
    }
}

# -------------------- Simulation --------------------
simulate_min <- function(rho = 0, g = 1, times = 0:500, params = list(), init = list()) {
    p <- modifyList(list(
        b=6.85e-5, delta=3.8e-5,
        omega=1e-4, gamma=0.4, delta_I=0.05,
        beta=1.2e-7, h_V=1e-5,
        sigma=50, delta_V=0.2,
        mu=8e-5, tau=80, delta_P=0.20,
        lambda=0.2,
        epsilon=0.1,
        theta=1, Z0=1,
        rho=rho, g=g,
        kappa=log(2)/(12/24)
    ), params)

    p <- unlist(p)   # <-- key step: force numerics

    y0 <- modifyList(list(
        S=1e6, I=0, R=0,
        VF=1, VZ=0, VI=0, P=4e3
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
times = 0:100
params_update <- list(b=0, omega = 0, mu = 2e-5, gamma = 0.6, delta_V = 0.3, epsilon = 0.05)
#params_update <- list(b=0, mu=6e-5, omega = 0, delta_V = 1e-1, c_lambda = 10, gamma = 0.4)
sim_no  <- simulate_min(rho=0.5, g=0, times=times, params = params_update)
sim_yes <- simulate_min(rho=0.5, g=0.9, times=times, params = params_update)

bind_rows(
    expected_daily(sim_no)  %>% mutate(model="No protection"),
    expected_daily(sim_yes) %>% mutate(model="Protection")
) %>%
    ggplot(aes(time, daily, color=model)) +
    geom_line(linewidth=0.8) +
    scale_y_continuous(labels=scales::comma) +
    labs(x="Time (days)", y="Expected daily infections",
         title="Daily infections with vs without plankton-conferred protection") +
    theme_bw() + theme(legend.position="top")

bind_rows(
    sim_no  %>% mutate(model="No protection"),
    sim_yes %>% mutate(model="Protection")
) %>%
    #select(time, VF, VZ, VI, P, model) %>%
    pivot_longer(-c(model, time)) %>%
    mutate(name = factor(name, c("S", "I", "R", "VF", "VI", "VZ", "P"))) %>%
    ggplot() +
    geom_line(aes(time, value, color=model)) +
    facet_wrap(~name, scales = "free_y") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()


epsilons <- c(0.5, 0.1, 0.05, 0.01)
sims <- map(epsilons, ~simulate_min(rho=0.5, g=0.8, params=list(epsilon=.x)))
names(sims) <- epsilons

df <- bind_rows(lapply(names(sims), function(e) {
    expected_daily(sims[[e]]) %>% mutate(epsilon=e)
}))

ggplot(df, aes(time, daily, color=epsilon)) +
    geom_line() +
    labs(y="Daily infections", x="Time",
         title="Effect of zooplankton detachment rate (ε) on epidemic dynamics")
