#'

library(tidyverse)
library(deSolve)
library(ggh4x)
library(cowplot)
source(here::here("tools.R"))

# Define the model as a function ----
model <- function(t, state, parameters) {

    with(as.list(c(state, parameters)), {
        # Equations
        dS <- b * N - delta * S - beta * V * S + omega * R
        dI <- beta * V * S - (gamma + delta_I + delta) * I
        dR <- gamma * I - (omega + delta) * R
        dV <- sigma * I + rho * Z_V - c * lambda * V * Z_F - mu * V * P - delta_V * V
        dP <- tau * mu * V * P - delta_P * P
        dZ_F <- eta * alpha * X * (Z_F + Z_V) - lambda * V * Z_F - delta_Z * Z_F
        dZ_V <- lambda * V * Z_F - delta_Z * Z_V
        dX <- r * X * (1 - X / K) - alpha * X * (Z_V + Z_F)

        list(c(dS, dI, dR, dV, dP, dZ_F, dZ_V, dX))
    })
}

# Initial state values (example, replace with your initial conditions) ----
state0 <- c(
    S = 1000,
    I = 0,
    R = 0,
    V = 0,
    P = 0,
    Z_F = 0,
    Z_V = 0,
    X = 0
)

parameters0 <- list(
    b = 6.85e-5,
    N = 1000,
    delta = 3.8e-5,
    beta = 0.214,
    omega = 0.00092,
    gamma = 0.1,
    delta_I = 0.013,
    sigma = 10,
    rho = 100,
    mu = 1.4e-9,
    delta_V = 0.33,
    tau = 80,
    delta_P = .5,
    lambda = 0.005,
    c = 5e-7,
    eta = 0.6,
    alpha = 0.4,
    delta_Z = 0.06,
    r = 0.6,
    K = 1
)

times <- seq(0, 100, by = 0.1)

# Update parameters ----
# Parameter
parameters <- modifyList(
    parameters0,
    list(
        # Turn off human demography
        b = 0,
        delta = 0,
        omega = 0,
        gamma = 0.2,
        # V to I, infection
        beta = 0.2,
        # I to V, shedding
        sigma = 0.01,
        # V to Z_V, infection
        lambda = 0.2,
        c = 1,
        # Z_V to V, shedding
        rho = 3,
        # V to P, predation
        mu = 1e-8,
        tau = 80
        #delta_P = 0
    )
)


## Init state
## Set X and P at equilibrium
X_star = parameters[["delta_Z"]] / (parameters[["eta"]] * parameters[["alpha"]])
Z_star = parameters[["r"]] / parameters[["alpha"]] * (1 - X_star / parameters[["K"]])

updates <- c(
    S = 1000,
    I = 0,
    R = 0,
    V = .1,
    P = 1e8,
    Z_F = Z_star,
    Z_V = 0,
    X = X_star
)
state <- state0
state[names(updates)] <- updates

# Run the simulation
out <- ode(y = state, times = times, func = model, parms = parameters) %>%
    clean_ode() %>%
    left_join(tb_types)

# Plot
p <- plot_one_out(out)
ggsave(here::here("plots/main-00.png"), p, width = 8, height = 4)

## Scenarios ----
parameters <- modifyList(parameters0, list(
    # Turn off human demography
    b = 0, delta = 0, omega = 0,
    #
    beta = 0.01, gamma = 0.2, c = 1,
    # Shedding
    sigma = .1, rho = 10,
    # Plankton
    alpha = 1e-5, eta = .5
))
parameters[["K"]] = 2*parameters[["delta_Z"]]/ ( parameters[["eta"]] * parameters[["alpha"]])

state <- state0
X_star = parameters[["delta_Z"]] / (parameters[["eta"]] * parameters[["alpha"]])
X_star
Z_star = parameters[["r"]] / parameters[["alpha"]] * (1 - X_star / parameters[["K"]])
Z_star
state["S"] <- 1e4; state["V"] <- 1; state["Z_F"] <- Z_star; state["X"] <- X_star


tb_mods <- tibble(
    scenario = 1:4,
    pars = list(
        modifyList(parameters,list(lambda = 0,mu = 0)),
        modifyList(parameters,list(lambda = 0,mu = 1e-6)),
        modifyList(parameters,list(lambda = 0.005,mu = 0)),
        modifyList(parameters,list(lambda = 0.005,mu = 1e-6))
    ),
    init = list(
        {updates <- c(P = 0); state[names(updates)] <- updates; state},
        {updates <- c(P = 1e9); state[names(updates)] <- updates; state},
        {updates <- c(P = 0); state[names(updates)] <- updates; state},
        {updates <- c(P = 1e9); state[names(updates)] <- updates; state}
    )
)

tb_results <- tb_mods %>%
    mutate(
        out = map2(pars, init, ~ode(y = .y, times = times, func = model, parms = .x)),
        tb = map(out, ~clean_ode(.x) %>% left_join(tb_types)),
        tb = map(tb, ~filter(.x, time <= 20))
    )

# Plot
p <- plot_grid(
    plot_one_out(tb_results$tb[[1]], "Scenario 1: SIRB"),
    plot_one_out(tb_results$tb[[2]], "Scenario 2: Phage only"),
    plot_one_out(tb_results$tb[[3]], "Scenario 3: Planktons only"),
    plot_one_out(tb_results$tb[[4]], "Scenario 4: Full model"),
    ncol = 1
)
p
ggsave(here::here("plots/mains.png"), p, width = 8, height = 12)





p <- plot_one_out(tb_results$tb[[1]], "Scenario 1: SIRB")
ggsave(here::here("plots/main-01.png"), p, width = 8, height = 4)
p <- plot_one_out(tb_results$tb[[2]], "Scenario 2: Phage")
ggsave(here::here("plots/main-02.png"), p, width = 8, height = 4)
p <- plot_one_out(tb_results$tb[[3]], "Scenario 3: Plankton")
ggsave(here::here("plots/main-03.png"), p, width = 8, height = 4)
p <- plot_one_out(tb_results$tb[[4]], "Scenario 4: Full model")
ggsave(here::here("plots/main-04.png"), p, width = 8, height = 4)

