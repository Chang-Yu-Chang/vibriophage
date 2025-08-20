#' Maity 2025

library(tidyverse)
library(deSolve)

clean_ode <- function (out) {
    as_tibble(out) %>%
        pivot_longer(-time) %>%
        mutate(time = as.numeric(time), value = as.numeric(value))
}

tb <- tibble(
    name = c("S", "I", "R", "B", "Z_B", "Z_F", "P"),
    type = c(rep("human", 3), rep("microbe", 4))
)


# Define the model
plankton_model <- function(time, state, parameters) {
    # Unpack state variables
    S <- state[["S"]]
    I <- state[["I"]]
    R <- state[["R"]]
    B <- state[["B"]]
    Z_B <- state[["Z_B"]]
    Z_F <- state[["Z_F"]]
    P <- state[["P"]]

    # Calculate derivatives
    with(as.list(parameters), {

        dS <- LAMBDA * N0 - beta * S * B / (h_b + B) - beta_z * S * Z_B / (h_z + Z_B) - mu * S + omega * R
        dI <- beta * S * B / (h_b + B) + beta_z * S * Z_B / (h_z + Z_B) - (gamma + mu + delta) * I
        dR <- gamma * I - (mu + omega) * R
        dB <- xi * I - d_b * B - c * sigma * B * Z_F / (h_m + B)
        dZ_B <- sigma * B * Z_F / (h_m + B) - d_z * Z_B
        dZ_F <- eta * alpha * P * (Z_F + Z_B) / (h_p + P) - d_z * Z_F - sigma * B * Z_F / (h_m + B)
        dP <- r_p * P * (1 - P/K) - alpha * P * (Z_F + Z_B) / (h_p + P)

        return(list(c(dS, dI, dR, dB, dZ_B, dZ_F, dP)))
    })
}

# Set parameters ----
parameters <- c(
    N0 = 0,              # Initial human population
    LAMBDA = 6.85e-5,    # Constant recruitment rate of human population. person/day
    beta = 0.214,        # Transmission rate via bacteria. 1/day
    beta_z = 1,          # Transmission rate via zooplankton. 1/day
    h_b = 1e9,           # Half saturation constant of bacterial transmission. cells/L
    h_z = 20,            # Half saturation constant of transmission via zooplankton. mg/L
    mu = 3.8e-5,         # Natural death rate of human. 1/day
    omega = 9.2e-4,      # Rate of immunity loss of recovered individuals. 1/day
    #omega = 0,
    delta = 0.013,       # Disease induced mortality of human. 1/day
    d_b = 0.33,          # Removal rate of the bacteria. 1/day
    d_z = 0.06,          # Death rate of zooplankton. 1/day
    gamma = 0.2,         # Recovery rate of infected human. 1/day
    xi = 10,             # Bacteria shedding rate of infected humans. cells/day/L per person
    sigma = 0.005,       # Rate of bacteria-zooplankton association. 1/day
    c = 5e7,             # Colonization coefficient of bacteria. cells/mg
    h_m = 1e7,           # Half saturation constant for bacteria-zooplankton association. cells/L
    eta = 0.6,           # Conversion coefficient
    alpha = 0.4,         # Predation rate. 1/day
    h_p = 0.6,           # Half saturation constant of phytoplankton growth. mg/L
    r_p = 0.5,           # Growth rate of phytoplankton. 1/day
    K = 1                # Carrying capacity of phytoplankton. mg/L
    # m = 0.3,            # Bacterial growth rate. Day^-1
    # Kv = 2.5e6,         # Bacterial carrying capacity. Cells per litter
    # beta = 100,         # Phage burst size. Virions per cell
    # gamma = 1.4e-9,     # Phage adsorption rate. Liters per virion per day
    # omega = 0.525,      # Phage decay rate. Virions per day
    # k = 4e7,            # Bacterial 50% infectious dose. Cells per day
    # l = 2.1e7,          # Phage 50% infectious dose. Virions per day
    # pi = 0.1,           # Severe disease rate. Day^-1
    # mun = 0.1,          # Phage-negative recovery rate. Day^-1
    # mup = 0.1,          # Phage-postive recovery rate. Day^-1
    # c = 10,             # Mean bacterial shed rate. Cells per liter per day
    # alpha = 1,          # Mean phage shed rate. Virions per cell
    # a = 7,              # Threshold parameter
    # delta = 1e-4        # Death/immigration rate. Day^-1
)
# Default parameter such that phi = Kv * beta * gamma / omega = 2/3

# Initial state ----
initial_state <- c(
    S = 1e3,      # Susceptible
    I = 0,        # Infected
    R = 0,        # Recovered
    B = 1e3,        # Free living bacteria,
    Z_B = 1,      # Bacteria-associated zooplankton
    Z_F = 1,      # Uncolonized zooplankton
    P = .5        # Phytoplankton
    # N = 1e6,        # Total human population
    # In = 0,         # Phage-negative infecteds
    # Ip = 0,         # Phage-positive infecteds
    # R = 0,          # Recovereds/deceased
    # V = 5*2.5*1e6,  # Vibrio density
    # P = 1e8         # Phage density
)

time <- seq(0, 500, by = 1)

# ----

ode(y = initial_state, times = time, func = plankton_model, parms = parameters) %>%
    clean_ode() %>%
    left_join(tb) %>%
    mutate(name = factor(name, c("S", "I", "R", "B", "Z_B", "Z_F", "P"))) %>%
    ggplot() +
    geom_line(aes(x = time, y = value, color = name)) +
    facet_wrap(~name, scales = "free_y") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()







