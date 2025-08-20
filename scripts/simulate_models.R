# simulate_models.R
# Models:
# M0 = SIRB (SIR + free-living Vibrio VF only)
# M1 = Plankton only (VF, VZ, ZF, ZV, X; no phage)
# M2 = Phage only (VF, P; no plankton, no VZ)
# M3 = Both (VF, VZ, P, ZF, ZV, X) with ONE attached pool VZ (no k-stage protection)

suppressPackageStartupMessages({
    library(deSolve)
    library(tidyverse)
})

# -------- Scenario-specific initial states --------
# You can still override any field with simulate_model(..., init = list(...))
scenario_init <- function(model) {
    base <- list(
        # Population
        S = 9.99e4, I = 5, R = 0,
        # Environment
        VF = 50, VZ = 0, P = 0,
        ZF = 0, ZV = 0, X = 0
    )

    add <- switch(
        model,
        # M0: SIRB only (free Vibrio present, no phage/plankton)
        "M0" = list(S=9.99e5, I=3, R=0, VF=50,  VZ=0, P=0,    ZF=0,  ZV=0,   X=0),

        # M1: Plankton-only (plankton bloom present; some VZ already attached)
        "M1" = list(S=9.99e5, I=5, R=0, VF=80,  VZ=3, P=0,    ZF=6,  ZV=0.3, X=1.5),

        # M2: Phage-only (environmental phage present; no plankton, VZ = 0)
        "M2" = list(S=9.99e5, I=5, R=0, VF=80,  VZ=0, P=4e5,  ZF=0,  ZV=0,   X=0),

        # M3: Both phage + plankton (moderate bloom + detectable phage)
        "M3" = list(S=9.99e5, I=5, R=0, VF=80,  VZ=3, P=4e5,  ZF=6,  ZV=0.3, X=1.5)
    )

    modifyList(base, add)
}

# ---------------- RHS factory ----------------
# We always return derivatives in this order:
# y = c(S, I, R, VF, VZ, P, ZF, ZV, X)
rhs_factory <- function(model = c("M0","M1","M2","M3")) {
    model <- match.arg(model)

    function(t, y, p) {
        with(as.list(c(y, p)), {
            N <- S + I + R

            # ------- Force of infection (saturating) -------
            # Veff = VF (+ theta * VZ when VZ is present)
            #theta_eff <- if (model %in% c("M1","M3")) theta else 0
            Veff <- VF + theta * VZ #* ifelse(model %in% c("M1","M3"), VZ, 0)
            lambda_H <- beta * Veff / (1 + h_V * Veff)

            dS <- b*N + omega*R - delta*S - lambda_H * S
            dI <- lambda_H * S - (gamma + delta_I + delta)*I
            dR <- gamma*I - (omega + delta)*R

            # ------- Defaults for ecological pieces -------
            dVF <- dVZ <- dP <- dZF <- dZV <- dX <- 0

            # ------- Vibrio (VF) base dynamics (all models) -------
            # SIRB core: shed from I, washout
            dVF <- dVF + sigma*I - delta_V*VF

            # ------- Plankton coupling (M1, M3) -------
            if (model %in% c("M1","M3")) {
                # Zooplankton & phyto
                dZF <- eta*alpha*X*(ZF + ZV) - lambda*VF*ZF + kappa_Z*ZV - delta_Z*ZF
                dZV <- lambda*VF*ZF - kappa_Z*ZV - delta_Z*ZV
                dX  <- r*X*(1 - X/K) - alpha*X*(ZF + ZV)

                # Colonization VF -> VZ and detachment VZ -> VF
                col_to_Z <- c_lambda*lambda*VF*ZF
                dVF <- dVF - col_to_Z + epsilon*VZ
                dVZ <- dVZ + col_to_Z - epsilon*VZ - delta_V*VZ
            }

            # ------- Phage dynamics (M2, M3) -------
            if (model %in% c("M2","M3")) {
                # Loss of VF to phage; in M3, attached VZ can also be hit with susceptibility g
                phage_loss_VF <- mu*VF*P
                phage_loss_VZ <- if (model == "M3") g*mu*VZ*P else 0
                dVF <- dVF - phage_loss_VF
                dVZ <- dVZ - phage_loss_VZ

                dP <- tau * (phage_loss_VF + phage_loss_VZ) - delta_P*P
            }

            # Return in fixed order matching the state vector
            list(c(dS, dI, dR, dVF, dVZ, dP, dZF, dZV, dX))
        })
    }
}

# --------------- Simulation wrapper ---------------
simulate_model <- function(model = c("M0","M1","M2","M3"),
                           times = 0:180,
                           params = list(),
                           init   = list()) {
    model <- match.arg(model)
    rhs   <- rhs_factory(model)

    # Default params (stable numerics; tune/override for data)
    p <- modifyList(list(
        # Human
        b=6.85e-5, delta=3.8e-5, beta=7e-8, h_V=1e-5,
        omega=9.2e-4, gamma=0.2, delta_I=0.013,
        # Vibrio & phage
        sigma=50, delta_V=0.28,
        mu=4e-8, tau=80, delta_P=0.15,
        # Plankton
        lambda=0.25, c_lambda=1.5, epsilon=0.5,
        eta=0.6, alpha=0.4, delta_Z=0.06, r=0.5, K=1,
        # Infectivity multipliers
        theta=0.7, g=0.6,
        # Zooplankton transition
        kappa_Z=12
    ), params)

    # Default initials (superset, unused compartments can be 0)
    y0 <- modifyList(scenario_init(model), init)

    # Single fixed state vector layout
    y <- c(S=y0$S, I=y0$I, R=y0$R, VF=y0$VF, VZ=y0$VZ, P=y0$P, ZF=y0$ZF, ZV=y0$ZV, X=y0$X)

    out <- ode(y = y, times = times, func = rhs, parms = p, method = "lsoda") %>% as_tibble()
    out$model <- model
    attr(out, "params") <- p
    out
}
