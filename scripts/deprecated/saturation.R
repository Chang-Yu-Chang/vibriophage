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
        dS <- b * N - delta * S - beta * V * S / (h_b + V) + omega * R
        dI <- beta * V * S / (h_b + V) - (gamma + delta_I + delta) * I
        dR <- gamma * I - (omega + delta) * R
        dV <- sigma * I + rho * Z_V - (c * lambda * V * Z_F) / (h_m + V) - (mu * V * P) / (h_P + P) - delta_V * V
        dP <- tau * mu * V * P - delta_P * P
        dZ_F <- (eta * alpha * X * (Z_F + Z_V)) / (h_X + X) - (lambda * V * Z_F) / (h_m + V) - delta_Z * Z_F
        dZ_V <- (lambda * V * Z_F) / (h_m + V) - delta_Z * Z_V
        dX <- r * X * (1 - X / K) - (alpha * X * (Z_V + Z_F)) / (h_X + X)


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
    gamma = 0.2,
    delta_I = 0.013,
    sigma = 1e-4,
    rho = 1e-2,
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
    K = 1,
    h_b = 100,
    h_P = 1,
    h_m = 1e3,
    h_X = 0.6
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
        lambda = 0.1,
        c = 1,
        # Z_V to V, shedding
        rho = 3,
        # V to P, predation
        mu = 0,
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
    V = 1e3,
    P = 1e3,
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
p

ggsave(here::here("plots/main-00.png"), p, width = 8, height = 4)

