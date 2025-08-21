# Minimal SIR–Vibrio–Phage–Zooplankton Model

This model couples human SIR dynamics with *Vibrio cholerae*, phage, and a zooplankton refuge. It formalizes two ideas:

1. **Refuge:** *Vibrio* attached to zooplankton (VZ) are shielded from phage;
2. **Temporary protection after detachment:** a fraction of released cells enter a transient, free-living **immune** compartment (VI) before reverting to ordinary free-living *Vibrio* (VF).

The key control is **$g\in[0,1]$**, which sets **how immune VI is to phage**:

* $g=1$: VI fully immune (no phage killing of VI)
* $g=0$: VI not immune (VI killed like VF)

We typically initialize with **zooplankton-associated *Vibrio*** (VZ) and **resident phage** present; VZ do not directly infect humans in the minimal presentation (set $\theta=0$, but this can be relaxed).

---

## State variables

* $S,I,R$: susceptible, infected, recovered humans
* $V_F$: free-living *Vibrio*
* $V_Z$: zooplankton-associated *Vibrio* (refuge)
* $V_I$: released, **temporarily immune** free-living *Vibrio*
* $P$: lytic phage

---

## Force of infection (FOI)

Effective *Vibrio* for infection:

$$
V_{\text{eff}} \;=\; V_F \;+\; V_I \;+\; \theta\,V_Z,
$$

and (minimal linear form used in code)

$$
\lambda_H \;=\; \beta\, V_{\text{eff}}.
$$

*(If desired, this can be replaced by a saturating or type-III form; see commented lines in the code.)*

---

## Model equations

Human SIR:

$$
\begin{aligned}
\frac{dS}{dt} &= b\,N + \omega R - \delta S - \lambda_H S, \\
\frac{dI}{dt} &= \lambda_H S - (\gamma + \delta_I + \delta) I, \\
\frac{dR}{dt} &= \gamma I - (\omega + \delta) R,\qquad N=S+I+R.
\end{aligned}
$$

Vibrio on zooplankton (refuge; **no** phage killing here):

$$
\frac{dV_Z}{dt} \;=\; \lambda\,V_F\,Z_0 \;-\; \epsilon\,V_Z \;-\; \delta_V\,V_Z.
$$

Released, **temporarily immune** Vibrio (phage loss scaled by $1-g$):

$$
\frac{dV_I}{dt} \;=\; \epsilon\,V_Z \;-\; \kappa\,V_I \;-\; (1-g)\,\mu\,V_I\,P \;-\; \delta_V\,V_I.
$$

Free-living Vibrio (fully phage-susceptible):

$$
\frac{dV_F}{dt} \;=\; \sigma\,I \;+\; \kappa\,V_I \;-\; \lambda\,V_F\,Z_0 \;-\; \mu\,V_F\,P \;-\; \delta_V\,V_F.
$$

Phage (replicate on $V_F$ and on $V_I$ only via its **residual** susceptibility $1-g$; no growth on $V_Z$):

$$
\frac{dP}{dt} \;=\; \tau\,\mu\left( V_F + (1-g)\,V_I \right) P \;-\; \delta_P\,P.
$$

---

## Parameters (subset; units per code)

* $\beta$ (L$^{-1}$,day$^{-1}$): infection scaling; $\theta$ (–): contribution of $V_Z$ to FOI (set to 0 in minimal runs).
* $\lambda$ (day$^{-1}$), $Z_0$ (–): VF–zooplankton encounters with constant zooplankton availability.
* $\epsilon$ (day$^{-1}$): **detachment** from $V_Z$ into $V_I$.
* $\kappa$ (day$^{-1}$): **loss of protection** $V_I \rightarrow V_F$ (mean protected time $\approx 1/\kappa$).
* $g \in [0,1]$: **immunity strength of $V_I$** against phage ( $g=1$ immune; $g=0$ none ).
* $\mu$ (L virion$^{-1}$ day$^{-1}$): phage adsorption; $\tau$: burst size; $\delta_P$ (day$^{-1}$): phage decay.
* $\sigma$ (cells host$^{-1}$ day$^{-1}$): shedding from infected humans; $\delta_V$ (day$^{-1}$): *Vibrio* loss/washout.
* $b,\delta,\omega,\gamma,\delta_I$: human demography and recovery.

---

## How we simulate (one paragraph)

We implement the ODEs in R with **deSolve**, using a constant zooplankton availability $Z_0$ to keep the ecology minimal. Runs typically start with a **zooplankton-associated seed** ($V_Z>0$) and **resident phage** present; humans start largely susceptible. We compare two scenarios by adjusting a **single parameter**: **$g=0$** (no post-detachment immunity; $V_I$ phage-susceptible) versus **$g=1$** (full immunity; $V_I$ phage-immune). The plotting helper `expected_daily()` uses the **same FOI** as the ODE to compute expected daily infections $\lambda_H S$.

---

## Main qualitative results (what the model shows)

* **No immunity ($g=0$)**: With resident phage, *all* free-living *Vibrio* (VF and VI) are phage prey. Phage expand and purge *Vibrio*, so **infections crash to zero** after a small transient—**no recurring outbreaks**.
* **With immunity ($g=1$)**: VI forms a **transient free-living refuge** that reseeds VF ($V_I \xrightarrow{\kappa} V_F$). Phage bloom on VF and crash it, but VI persists (protected) and **restarts the cycle**, yielding **predator–prey oscillations** (VF–P) and **recurring infection peaks**.
* **Role of detachment ($\epsilon$)**: Smaller $\epsilon$ (slower release from $V_Z$) acts like a slow-release reservoir, **broadening conditions** for multiple waves; larger $\epsilon$ shortens protection and favors a single-peak crash.

> **Take-home:** In a phage-present environment, **zooplankton refuge alone is not enough** to generate recurrent outbreaks if VZ does not directly infect; **post-detachment immunity (VI, controlled by $g$) is necessary** to sustain epidemic cycles.

---
