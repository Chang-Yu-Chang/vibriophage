#' Jenson 2006

library(tidyverse)
library(deSolve)

# Define the model
phage_model <- function(time, state, parameters) {
    # Unpack state variables
    N <- state[1]
    In <- state[2]
    Ip <- state[3]
    R <- state[4]
    V <- state[5]
    P <- state[6]
    S = N - In - Ip - R

    # Unpack parameters
    m = parameters[1]
    Kv = parameters[2]
    beta = parameters[3]
    gamma = parameters[4]
    omega = parameters[5]
    k = parameters[6]
    l = parameters[7]
    pi = parameters[8]
    mun = parameters[9]
    mup = parameters[10]
    c = parameters[11]
    alpha = parameters[12]
    a = parameters[13]
    delta = parameters[14]

    # Calculate derivatives
    # dprey <- alpha * prey - beta * prey * predator
    # dpredator <- delta * prey * predator - gamma * predator

    Ca <- function (a) 2^(1/a) - 1

    dS <- -pi * S * (V / (Ca(a)*k + V))^a - delta * S + delta * N
    dIn <- pi * S * (l / (l + P)) * (V / (Ca(a)*k + V))^a - (mun + delta) * In
    dIp <- pi * S * (P / (l + P)) * (V / (Ca(a)*k + V))^a - (mup + delta) * Ip
    dR <- mun * In + mup * Ip - delta * R
    dV <- V * (m * (1 - V/Kv) - gamma * P) + c * (In + Ip)
    dP <- P * (beta * gamma * V - omega) + alpha * c * Ip


    return(list(c(dS, dIn, dIp, dR, dV, dP)))
}

# Set parameters ----
parameters <- c(
    m = 0.3,            # Bacterial growth rate. Day^-1
    Kv = 2.5e6,         # Bacterial carrying capacity. Cells per litter
    beta = 100,         # Phage burst size. Virions per cell
    gamma = 1.4e-9,     # Phage adsorption rate. Liters per virion per day
    omega = 0.525,      # Phage decay rate. Virions per day
    k = 4e7,            # Bacterial 50% infectious dose. Cells per day
    l = 2.1e7,          # Phage 50% infectious dose. Virions per day
    pi = 0.1,           # Severe disease rate. Day^-1
    mun = 0.1,          # Phage-negative recovery rate. Day^-1
    mup = 0.1,          # Phage-postive recovery rate. Day^-1
    c = 10,             # Mean bacterial shed rate. Cells per liter per day
    alpha = 1,          # Mean phage shed rate. Virions per cell
    a = 7,              # Threshold parameter
    delta = 1e-4        # Death/immigration rate. Day^-1
)
# Default parameter such that phi = Kv * beta * gamma / omega = 2/3

# Initial state ----
initial_state <- c(
    N = 1e6,        # Total human population
    In = 0,         # Phage-negative infecteds
    Ip = 0,         # Phage-positive infecteds
    R = 0,          # Recovereds/deceased
    V = 5*2.5*1e6,  # Vibrio density
    P = 1e8         # Phage density
)

# Time vector for simulation ----
time <- seq(0, 100, by = 1)

# Run the simulation ----
clean_ode <- function (out) {
    as_tibble(out) %>%
    pivot_longer(-time) %>%
    mutate(time = as.numeric(time), value = as.numeric(value))
}

tb_mods <- tibble(
    hypo = c(rep("bloom", 2), rep("instability", 4)),
    mod = c("bloom P+", "bloom P-", "instability 1", "instability 2", "instability 3", "instability 4"),
    pars = list(
        parameters, parameters,
        {parameters["Kv"] = 2.5e6 * 9/4; parameters["omega"] = 0.525; parameters["a"] = 7; parameters}, # phi = 1.5
        {parameters["Kv"] = 2.5e6 * 9/4; parameters["omega"] = 0.394; parameters["a"] = 7; parameters}, # phi = 2
        {parameters["Kv"] = 2.5e6 * 9/4; parameters["omega"] = 0.263; parameters["a"] = 5; parameters}, # phi = 3
        {parameters["Kv"] = 2.5e6 * 9/4; parameters["omega"] = 0.158; parameters["a"] = 3; parameters} # phi = 5
    ),
    init = list(
        initial_state,
        {initial_state["P"] = 0; initial_state},
        {initial_state["P"] = 1; initial_state},
        {initial_state["P"] = 1; initial_state},
        {initial_state["P"] = 1; initial_state},
        {initial_state["P"] = 1; initial_state}
    )
) %>%
    mutate(
        out = map2(pars, init, ~ode(y = .y, times = time, func = phage_model, parms = .x)),
        tb = map(out, clean_ode)
    )

# Plot ----
p <- tb_mods %>%
    select(hypo, mod, tb) %>%
    unnest(tb) %>%
    filter(name %in% c("In", "Ip")) %>%
    #filter(name %in% c("Ip")) %>%
    group_by(hypo, mod, time) %>%
    summarize(value = sum(value)) %>%
    ggplot() +
    geom_line(aes(x = time, y = value, color = mod), linewidth = 1, alpha = .8) +
    facet_wrap(~hypo, scales = "free_y") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(here::here("plots/Jensen2006.png"), p, width = 10, height = 5)



if (F) {

tb <- ode(y = initial_state, times = time, func = phage_model, parms = parameters) %>%
    as_tibble %>%
    pivot_longer(-time) %>%
    mutate(time = as.numeric(time), value = as.numeric(value))

tb %>%
    #filter(name %in% c("V", "P"))
    filter(name %in% c("In", "Ip")) %>%
    group_by(time) %>%
    summarize(value = sum(value)) %>%
    ggplot() +
    geom_line(aes(x = time, y = value)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

}
