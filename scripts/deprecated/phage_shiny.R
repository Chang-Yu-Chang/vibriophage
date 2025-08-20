#'

library(tidyverse)
library(shiny)
library(deSolve)
library(cowplot)

# Define the model
phage_model <- function(time, state, parameters) {
    S <- state[1]
    In <- state[2]
    Ip <- state[3]
    R <- state[4]
    V <- state[5]
    P <- state[6]
    #S = N - In - Ip - R

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
    N = parameters[15]

    Ca <- function (a) 2^(1/a) - 1

    dS <- -pi * S * (V / (Ca(a)*k + V))^a - delta * S + delta * N
    dIn <- pi * S * (l / (l + P)) * (V / (Ca(a)*k + V))^a - (mun + delta) * In
    dIp <- pi * S * (P / (l + P)) * (V / (Ca(a)*k + V))^a - (mup + delta) * Ip
    dR <- mun * In + mup * Ip - delta * R
    dV <- V * (m * (1 - V/Kv) - gamma * P) + c * (In + Ip)
    dP <- P * (beta * gamma * V - omega) + alpha * c * Ip

    return(list(c(dS, dIn, dIp, dR, dV, dP)))
}
ode_equations <- c(
    "$$
    \\begin{aligned}
    \\frac{dS}{dt} = -\\pi S \\left( \\frac{V}{2^{1/a}k + V} \\right)^a - \\delta S + \\delta N \\\\
    \\frac{dI_n}{dt} = \\pi S \\left( \\frac{l}{l + P} \\right) \\left( \\frac{V}{2^{1/a}k + V} \\right)^a - (\\mu_n + \\delta) I_n \\\\
    \\frac{dI_p}{dt} = \\pi S \\left( \\frac{P}{l + P} \\right) \\left( \\frac{V}{2^{1/a}k + V} \\right)^a - (\\mu_p + \\delta) I_p \\\\
    \\frac{dR}{dt} = \\mu_n I_n + \\mu_p I_p - \\delta R \\\\
    \\frac{dV}{dt} = V \\left( m \\left( 1 - \\frac{V}{K_v} \\right) - \\gamma P \\right) + c \\left( I_n + I_p \\right) \\\\
    \\frac{dP}{dt} = P \\left( \\beta \\gamma V - \\omega \\right) + \\alpha c I_p \\\\
    \\end{aligned}
    $$
    "
)

# UI
ui <- fluidPage(
    titlePanel("Phage Model"),

    sidebarLayout(
        sidebarPanel(
            h3("Model Parameters"),
            sliderInput("m", expression("Bacterial growth rate (m):"), min = 0, max = 1, value = 0.3, step = 0.01),
            sliderInput("Kv", expression("Carrying capacity (K[v]):"), min = 1e6, max = 1e7, value = 2.5e6, step = 1e5),
            sliderInput("omega", expression("Phage decay rate (ω):"), min = 0, max = 1, value = 0.5, step = 0.01),
            sliderInput("a", expression("Threshold parameter (α):"), min = 1, max = 10, value = 7, step = 1),

            hr(), # Horizontal line separator

            h3("Initial State"),
            sliderInput("V", expression("Initial Vibrio density (V), exponential scale:"), min = log10(1), max = log10(1e10), value = log10(5*2.5*1e6), step = .1),
            sliderInput("P", expression("Initial Phage density (P), exponential scale:"), min = log10(1), max = log10(1e10), value = log10(1e8), step = .1)),

        mainPanel(
            fluidRow(
                column(8, plotOutput("phagePlot", height = "600px")),  # Plot on the left
                column(4, uiOutput("equations"))  # Equations on the right
            )
        )
    )
)

# Server logic
server <- function(input, output) {
    output$phagePlot <- renderPlot({
        parameters <- c(
            m = input$m,         # Bacterial growth rate. Day^-1
            Kv = input$Kv,       # Bacterial carrying capacity. Cells per litter
            beta = 100,          # Phage burst size. Virions per cell
            gamma = 1.4e-9,      # Phage adsorption rate. Liters per virion per day
            omega = input$omega, # Phage decay rate. Virions per day
            k = 4e7,             # Bacterial 50% infectious dose. Cells per day
            l = 2.1e7,           # Phage 50% infectious dose. Virions per day
            pi = 0.1,            # Severe disease rate. Day^-1
            mun = 0.1,           # Phage-negative recovery rate. Day^-1
            mup = 0.1,           # Phage-postive recovery rate. Day^-1
            c = 10,              # Mean bacterial shed rate. Cells per liter per day
            alpha = 1,           # Mean phage shed rate. Virions per cell
            a = input$a,         # Threshold parameter
            delta = 1e-4,        # Death/immigration rate. Day^-1
            N = 1e6              # Total human population size
        )

        initial_state <- c(
            S = 1e6,        # Susceptible
            In = 1,         # Phage-negative infecteds
            Ip = 1,         # Phage-positive infecteds
            R = 0,          # Recovereds
            V = 10^input$V, # Vibrio density
            P = 10^input$P  # Phage density
        )

        # Define the time vector
        time <- seq(0, 100, by = 1)

        # Run the simulation
        out <- ode(y = initial_state, times = time, func = phage_model, parms = parameters)

        # Clean the output into a tidy format
        clean_ode <- function(out) {
            as_tibble(out) %>%
                pivot_longer(-time) %>%
                mutate(time = as.numeric(time), value = as.numeric(value))
        }

        tb <- clean_ode(out)

        # Create the plot
        p1 <- tb %>%
            filter(name %in% c("In", "Ip")) %>%
            group_by(time) %>%
            summarize(value = sum(value)) %>%
            ggplot() +
            geom_line(aes(x = time, y = value), linewidth = 1, alpha = 0.8) +
            theme_bw() +
            labs(title = "Number of infected (Ip + In)", x = "Time (days)", y = "Cases")
        # p2 <- tb %>%
        #     filter(name %in% c("S", "R")) %>%
        #     ggplot() +
        #     geom_line(aes(x = time, y = value, color = name), linewidth = 1) +
        #     coord_cartesian(clip = "off") +
        #     theme_bw() +
        #     labs(title = "Number of susceptible (S) and recovered (R)", x = "Time (days)", y = "Cases")
        p2 <- tb %>%
            filter(name %in% c("V", "P")) %>%
            ggplot() +
            geom_line(aes(x = time, y = value, color = name), linewidth = 1) +
            scale_color_manual(values = c(P = "steelblue", V = "darkred"), labels = c(P = "Phage", V = "Vibrio")) +
            scale_y_log10() +
            coord_cartesian(clip = "off") +
            theme_bw() +
            theme() +
            guides() +
            labs(
                title = "Vibrio and Phage density",
                x = "Time (days)",
                y = "Density (log)"
            )

        p <- plot_grid(p1, p2, align = "v", ncol = 1)
        print(p)
    })

    # Render equations
    output$equations <- renderUI({
            withMathJax(ode_equations[1])
    })
}

# Run the application
shinyApp(ui = ui, server = server)
