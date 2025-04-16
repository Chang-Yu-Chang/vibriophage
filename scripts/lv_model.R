#'

library(tidyverse)
library(deSolve)

# Define the Lotka-Volterra model
lv_model <- function(time, state, parameters) {
    # Unpack state variables
    prey <- state[1]
    predator <- state[2]

    # Unpack parameters
    alpha <- parameters[1]  # Prey growth rate
    beta <- parameters[2]   # Predator's consumption rate
    delta <- parameters[3]  # Predator reproduction rate
    gamma <- parameters[4]  # Predator death rate

    # Calculate derivatives
    dprey <- alpha * prey - beta * prey * predator
    dpredator <- delta * prey * predator - gamma * predator

    return(list(c(dprey, dpredator)))  # Return the rates of change
}

# Set parameters
parameters <- c(alpha = 0.1, beta = 0.02, delta = 0.01, gamma = 0.1)

# Initial state: number of prey and predators
initial_state <- c(prey = 40, predator = 9)

# Time vector for simulation
time <- seq(0, 200, by = 1)

# Run the simulation
output <- ode(y = initial_state, times = time, func = lv_model, parms = parameters)

# Convert output to a data frame for easier handling
output_df <- as.data.frame(output)

# Plot the results
ggplot(data = output_df, aes(x = time)) +
    geom_line(aes(y = prey, color = "Prey")) +
    geom_line(aes(y = predator, color = "Predator")) +
    labs(title = "Lotka-Volterra Model", x = "Time", y = "Population Size") +
    scale_color_manual(values = c("Prey" = "blue", "Predator" = "red")) +
    theme_minimal()
