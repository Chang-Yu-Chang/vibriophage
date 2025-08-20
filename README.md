# vibriophage — SIRB + (phage | plankton) minimal models

This repository simulates cholera dynamics by coupling a human **SIR** model to environmental **Vibrio** and, depending on the variant, to **vibriophages** and **zooplankton**. We use a short, single **zooplankton-associated Vibrio** pool \(V_Z\) (no k-stage chain); attached cells detach back to \(V_F\) at rate \(\epsilon\). The force of infection (FOI) for humans saturates with Vibrio concentration.

## State variables
- Humans: \(S, I, R\) with \(N=S+I+R\).
- Vibrio & phage: \(V_F\) (free-living Vibrio), \(V_Z\) (zooplankton-associated Vibrio), \(P\) (phage).
- Plankton: \(Z_F\) (free zooplankton), \(Z_V\) (Vibrio-associated zooplankton), \(X\) (phytoplankton).

## Model ladder
- **M0 (SIRB)**: SIR + \(V_F\) (no \(V_Z\), no \(P\)).
- **M1 (Plankton-only)**: adds \(V_Z, Z_F, Z_V, X\); no phage.
- **M2 (Phage-only)**: adds \(P\); no plankton (\(V_Z=0\)).
- **M3 (Both)**: includes \(V_Z\) **and** \(P\). \(V_Z\) can be partially phage-susceptible via parameter \(g\).

## Dynamics (M3 "both" model; specialized models drop the absent terms)
Human FOI (saturating):
\[
\lambda_H(t) \;=\; \frac{\beta \, (V_F + \theta V_Z)}{1 + h_V \, (V_F + \theta V_Z)}.
\]

Humans:
\[
\begin{aligned}
\dot S &= bN + \omega R - \delta S - \lambda_H S,\\
\dot I &= \lambda_H S - (\gamma + \delta_I + \delta) I,\\
\dot R &= \gamma I - (\omega + \delta) R.
\end{aligned}
\]

Vibrio & phage:
\[
\begin{aligned}
\dot V_F &= \sigma I \;-\; \delta_V V_F \;-\; \underbrace{\mu V_F P}_{\text{phage}} \;-\; \underbrace{c_\lambda \lambda V_F Z_F}_{\text{colonize}} \;+\; \underbrace{\epsilon V_Z}_{\text{detach}},\\
\dot V_Z &= c_\lambda \lambda V_F Z_F \;-\; \epsilon V_Z \;-\; \delta_V V_Z \;-\; \underbrace{g\,\mu V_Z P}_{\text{phage on attached (M3 only)}},\\
\dot P   &= \tau \mu  ( V_F + g\, V_Z ) P \;-\; \delta_P P.
\end{aligned}
\]

Zooplankton & phytoplankton:
\[
\begin{aligned}
\dot Z_F &= \eta \alpha X (Z_F + Z_V)  + \kappa_Z Z_V - ( \lambda V_F + \delta_Z) Z_F,\\
\dot Z_V &= \lambda V_F Z_F - ( \kappa_Z + \delta_Z) Z_V,\\
\dot X   &= r X \Bigl(1 - \frac{X}{K}\Bigr) - \alpha X (Z_F + Z_V).
\end{aligned}
\]

Special cases:

- **M0**: set \(Z_F=Z_V=X=P=V_Z=0\); \(\dot V_F = \sigma I - \delta_V V_F\).
- **M1**: set \(P=0\) and \(g=0\).
- **M2**: set \(Z_F=Z_V=X=V_Z=0\).
- **M3**: as above (both).

## Parameters (key ones)
- Human: \(\beta\) (scale), \(h_V\) (half-saturation), \(\gamma, \omega, b, \delta, \delta_I\).
- Vibrio: \(\sigma\) (shedding), \(\delta_V\) (loss).
- Phage: \(\mu\) (adsorption), \(\tau\) (burst), \(\delta_P\) (decay), \(g\in[0,1]\) (phage access to \(V_Z\); M3).
- Plankton: \(\lambda\) (VF–ZF encounter), \(c_\lambda\) (colonization efficiency), \(\epsilon\) (detachment), \(\eta,\alpha,\delta_Z,r,K,\kappa_Z\).
- Infectivity mixing: \(\theta\) (contribution of \(V_Z\) to FOI).

Units are kept consistent with day\(^{-1}\) rates; \(\beta\) and \(h_V\) together control the scale/saturation of human FOI.

## From model to data (what we compare)
We compare **daily incidence** (new infections), not \(I(t)\):
\[
E[Y_t] \;=\; \rho \;\lambda_H(t) S(t),
\]
where \(\rho\) is a reporting/scale factor fit by maximum likelihood for each model.

## Scoring
- **RMSE** on daily cases (lower is better).
- **Negative Binomial log score** with fixed dispersion \(\phi\) (higher is better).  
Both are implemented in `observe_and_score.R`:
```r
ex <- expected_daily_cases(sim)
rho <- fit_rho_daily(obs_daily, sim, start_date)
rmse_daily(obs$cases, rho * ex$mu)
nb_logscore_daily(obs$cases, rho * ex$mu, phi = 10)
