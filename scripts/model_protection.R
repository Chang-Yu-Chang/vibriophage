# figure2_minimal.R
# Two tiny models (No-Protection vs Protection) + B1 (time series) + B2 (robustness bars)

suppressPackageStartupMessages({
    library(deSolve)
    library(tidyverse)
    library(patchwork)
})

# -------------------- Core ODEs --------------------
# States:
#  - No-Protection model: y = c(S,I,R, VF, VZ, P)
#  - Protection model:    y = c(S,I,R, VF, VZ, VI, P)
#
# VZ = zooplankton-associated Vibrio (attached)
# VI = released, immune free-living Vibrio that loses protection (12h half-life)
#
# FOI: lambda_H = beta * Veff / (1 + h_V * Veff), Veff = VF + theta*(VZ + VI)
# Colonization uses *constant* zooplankton availability Z0 to keep things minimal.

rhs_factory <- function(protection = TRUE) {
    function(t, y, p) {
        with(as.list(c(y, p)), {
            N <- S + I + R
            VZx <- if (protection) VZ else VZ
            VIx <- if (protection) VI else 0
            Veff <- VF + theta * (VZx + VIx)
            lambda_H <- beta * Veff / (1 + h_V * Veff)

            dS <- b*N + omega*R - delta*S - lambda_H*S
            dI <- lambda_H*S - (gamma + delta_I + delta)*I
            dR <- gamma*I - (omega + delta)*R

            # Vibrio base & colonization
            col_to_Z <- c_lambda * lambda * VF * Z0  # Z0 is constant availability
            dVF <- sigma*I - delta_V*VF - col_to_Z
            dVZ <- col_to_Z - epsilon*VZ - delta_V*VZ

            # Detachment path differs by model:
            if (protection) {
                # VZ -> VI (immune), then VI -> VF at kappa_I
                dVI <- epsilon*VZ - kappa_I*VI - delta_V*VI         # immune: no phage loss
                dVF <- dVF + kappa_I*VI                              # return flow to VF
            } else {
                # No protection: detachment directly to VF
                dVI <- 0
                dVF <- dVF + epsilon*VZ
            }

            # Phage dynamics:
            # - No-Protection: phage kill both VF and VZ (VZ susceptibility = g=1)
            # - Protection:    phage kill VF only (VZ/VI immune during window)
            if (protection) {
                loss_VF <- mu*VF*P
                dVF <- dVF - loss_VF
                dP  <- tau*loss_VF - delta_P*P
            } else {
                loss_VF <- mu*VF*P
                loss_VZ <- mu*VZ*P
                dVF <- dVF - loss_VF
                dVZ <- dVZ - loss_VZ
                dP  <- tau*(loss_VF + loss_VZ) - delta_P*P
            }

            if (protection) {
                list(c(dS, dI, dR, dVF, dVZ, dVI, dP))
            } else {
                list(c(dS, dI, dR, dVF, dVZ, dP))
            }
        })
    }
}

simulate_min <- function(protection = TRUE, times = 0:120, params = list(), init = list()) {
    # sensible defaults (minimal & stable)
    p <- modifyList(list(
        b=6.85e-5, delta=3.8e-5, omega=9.2e-4, gamma=0.2, delta_I=0.013,
        beta=1.2e-7, h_V=1e-5,
        sigma=50, delta_V=0.30,
        mu=3e-6, tau=80, delta_P=0.20,
        lambda=0.18, c_lambda=1.2, epsilon=0.50,
        theta=0.8,
        Z0=6,                    # constant zooplankton availability (units arbitrary)
        kappa_I = log(2)/0.5     # ~12 h half-life => 0.5 day
    ), params)

    if (protection) {
        y0 <- modifyList(list(
            S=9.99e5, I=5, R=0,
            VF=80, VZ=3, VI=0, P=4e5
        ), init)
        y <- c(S=y0$S, I=y0$I, R=y0$R, VF=y0$VF, VZ=y0$VZ, VI=y0$VI, P=y0$P)
    } else {
        y0 <- modifyList(list(
            S=9.99e5, I=5, R=0,
            VF=80, VZ=3, P=4e5
        ), init)
        y <- c(S=y0$S, I=y0$I, R=y0$R, VF=y0$VF, VZ=y0$VZ, P=y0$P)
    }

    rhs <- rhs_factory(protection)
    out <- ode(y=y, times=times, func=rhs, parms=p, method="lsoda") %>% as_tibble()
    attr(out, "params") <- p
    out
}

# expected daily infections (what we plot/score)
expected_daily <- function(sim_df, protection = TRUE) {
    p <- attr(sim_df, "params")
    VZ <- if ("VZ" %in% names(sim_df)) sim_df$VZ else 0
    VI <- if (protection && "VI" %in% names(sim_df)) sim_df$VI else 0
    Veff <- sim_df$VF + p$theta * (VZ + VI)
    lambda_H <- p$beta * Veff / (1 + p$h_V * Veff)
    tibble(time = sim_df$time, daily = lambda_H * sim_df$S)
}

# -------------------- Panel B1: time series --------------------
times <- 0:120
sim_no_prot <- simulate_min(protection = FALSE, times = times)
sim_prot    <- simulate_min(protection = TRUE,  times = times, params = list(kappa_I = log(2)/(12/24)))

b1_df <- bind_rows(
    expected_daily(sim_no_prot, protection = FALSE) %>% mutate(model = "No protection") %>%
        mutate(across(-model, as.numeric)),
    expected_daily(sim_prot,    protection = TRUE)  %>% mutate(model = "12h protection") %>%
        mutate(across(-model, as.numeric))
)

p_B1 <- ggplot(b1_df, aes(time, daily, color = model)) +
    geom_line(size = 1) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Time (days)", y = "Expected daily infections",
         title = "Minimal model: daily infections",
         subtitle = "With vs without 12h plankton-induced protection") +
    theme_bw() +
    theme(legend.position = "top")
p_B1
# -------------------- Panel B2: robustness bars --------------------
# Features from a daily-infections time series df = tibble(time, daily)
curve_features <- function(df, drop_window_days = 21) {  # shorter window (3 weeks)
    dt    <- median(diff(df$time))
    kdrop <- max(1, round(drop_window_days / dt))
    imax  <- which.max(df$daily)
    Imax  <- df$daily[imax]
    tpeak <- df$time[imax]
    Iafter <- df$daily[pmin(nrow(df), imax + kdrop)]
    crash <- if (Imax > 0) (Imax - Iafter) / Imax else NA_real_
    tibble(Imax = Imax, tpeak = tpeak, crash = crash)
}

# More permissive "data-like" test for DAILY infections
# Peak 0.01%–2% of N per day, peak between day 5–120, and ≥50% drop within 21 days
is_data_like <- function(feat, N,
                         p_lo = 1e-4, p_hi = 2e-2,    # 0.01%–2% of N per day
                         t_lo = 5, t_hi = 120,        # peak timing window
                         crash_ge = 0.50,             # at least 50% drop
                         Imin_abs = 50) {             # avoid trivial peaks
    peak_prop <- feat$Imax / N
    ok_peak   <- (peak_prop >= p_lo) && (peak_prop <= p_hi)
    ok_time   <- (feat$tpeak >= t_lo) && (feat$tpeak <= t_hi)
    ok_crash  <- (feat$crash >= crash_ge)
    ok_mass   <- (feat$Imax >= Imin_abs)
    ok_peak && ok_time && ok_crash && ok_mass
}

# Scan (μ, ε) grids; now uses N from the underlying simulation
scan_robust <- function(protection, mu_vals, eps_vals, params = list(), init = list(), times = 0:120) {
    grid <- expand_grid(mu = mu_vals, epsilon = eps_vals)
    purrr::pmap_dfr(grid, function(mu, epsilon) {
        sim <- simulate_min(protection = protection,
                            times = times,
                            params = modifyList(params, list(mu = mu, epsilon = epsilon)),
                            init = init)
        # compute N from the trajectory (S+I+R)
        Npop <- max(sim$S + sim$I + sim$R, na.rm = TRUE)
        ed   <- expected_daily(sim, protection = protection)
        feat <- curve_features(ed, drop_window_days = 21)
        tibble(mu = mu, epsilon = epsilon, data_like = is_data_like(feat, N = Npop))
    })
}
mu_vals  <- 10^seq(-8, -5, length.out = 14)   # 1e-8 ... 1e-5
eps_vals <- seq(0.2, 1.2, length.out = 11)    # 0.2 ... 1.2 day^-1

res_no   <- scan_robust(FALSE, mu_vals, eps_vals)
res_yes  <- scan_robust(TRUE,  mu_vals, eps_vals)

sumtab <- bind_rows(
    res_no  %>% summarise(frac = mean(data_like)) %>% mutate(model = "No protection"),
    res_yes %>% summarise(frac = mean(data_like)) %>% mutate(model = "12h protection")
)

p_B2 <- ggplot(sumtab, aes(model, frac, fill = model)) +
    geom_col(width = 0.6) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("No protection"="#999999","12h protection"="#4CAF50")) +
    labs(x = NULL, y = "Share of (μ, ε) pairs that are data-like",
         title = "Robustness: protection widens the data-like region") +
    theme_bw() +
    theme(legend.position = "none")

# -------------------- (Optional) third panel: heatmap for protection case --------------------
res_yes$label <- ifelse(res_yes$data_like, "data-like", "other")
p_heat_opt <- ggplot(res_yes, aes(mu, epsilon, fill = label)) +
    geom_tile() +
    scale_x_log10() +
    scale_fill_manual(values = c("data-like"="#4CAF50", "other"="#f0f0f0")) +
    labs(x = expression(mu~"(adsorption)"), y = expression(epsilon~"(detach)"),
         title = "Parameter region (12h protection)") +
    theme_bw() + theme(legend.title = element_blank())

# -------------------- Assemble figure(s) --------------------
# Two-panel layout (B1 + B2):
p_two <- p_B1 / p_B2
print(p_two)

# Uncomment for the 3-panel version:
# p_three <- p_B1 / (p_B2 | p_heat_opt)
# print(p_three)
