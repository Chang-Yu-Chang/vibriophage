# vibriophage — SIRB + (phage | plankton) minimal models

This repository simulates cholera dynamics by coupling a human **SIR** model to environmental **Vibrio** and, depending on the variant, to **vibriophages** and **zooplankton**. We use a short, single **zooplankton-associated Vibrio** pool \(V_Z\) (no k-stage chain); attached cells detach back to \(V_F\) at rate \(\epsilon\). The force of infection (FOI) for humans saturates with Vibrio concentration.

Here’s your model description in clean **Markdown** format with equations in LaTeX blocks. You can paste this directly into your `README.md` and GitHub will render the math if you use a LaTeX renderer (e.g. with MathJax via GitHub Pages, Quarto, or JupyterBook).

---

# Model Equations

We track human hosts $(S,I,R)$, free Vibrio $(V_F)$, zooplankton-associated Vibrio $(V_Z)$, optionally immune Vibrio $(V_I)$, and phage $(P)$.

Total host population:

```math
N = S + I + R
```

**Force of infection (per susceptible host):**

```math
\lambda_H = \frac{\beta \, (V_F + \theta (V_Z + V_I))}{1 + h_V \, (V_F + \theta (V_Z + V_I))}
```

---

## Host dynamics

```math
\begin{aligned}
\frac{dS}{dt} &= bN + \omega R - \delta S - \lambda_H S \\
\frac{dI}{dt} &= \lambda_H S - (\gamma + \delta_I + \delta) I \\
\frac{dR}{dt} &= \gamma I - (\omega + \delta) R
\end{aligned}
```

---

## Vibrio dynamics (base terms)

Shedding into water:

```math
\sigma I
```

Colonization to zooplankton:

```math
\text{col\_to\_Z} = c_\lambda \, \lambda \, V_F Z_0
```

```math
\begin{aligned}
\frac{dV_F}{dt} &= \sigma I - \delta_V V_F - \text{col\_to\_Z} \\
\frac{dV_Z}{dt} &= \text{col\_to\_Z} - \epsilon V_Z - \delta_V V_Z
\end{aligned}
```

---

## Case 1. With zooplankton protection

When detaching, zooplankton-associated Vibrio enter a **protected state** $(V_I)$, immune to phage for $\sim 12h$:

```math
\begin{aligned}
\frac{dV_I}{dt} &= \epsilon V_Z - \kappa_I V_I - \delta_V V_I \\
\frac{dV_F}{dt} &+= \kappa_I V_I
\end{aligned}
```

Phage kill only free Vibrio:

```math
\begin{aligned}
\frac{dV_F}{dt} &-= \mu V_F P \\
\frac{dP}{dt} &= \tau \mu V_F P - \delta_P P
\end{aligned}
```

---

## Case 2. No protection

Detachment feeds directly back to free Vibrio:

```math
\frac{dV_F}{dt} += \epsilon V_Z
```

Phage kill both free and zooplankton-associated Vibrio:

```math
\begin{aligned}
\frac{dV_F}{dt} &-= \mu V_F P \\
\frac{dV_Z}{dt} &-= \mu V_Z P \\
\frac{dP}{dt} &= \tau \mu (V_F + V_Z) P - \delta_P P
\end{aligned}
```

---

⚠️ Notes:

* $V_I$ appears only in the **protection** scenario.
* In protection, phage kill **only $V_F$**.
* In no-protection, phage kill both $V_F$ and $V_Z$.
* Colonization ($\text{col\_to\_Z}$) always depends on $V_F$ and fixed $Z_0$.

---

👉 Do you also want me to add a **parameter table** (with biological interpretation and default values from your code) in the README? That would make it more self-contained for readers.
