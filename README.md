# vibriophage

This repository includes scripts to simulate the SIRB vibrio-phage-plankton dynamics. The main model is modified from two papers:

- [Jensen et al 2006](https://www.pnas.org/doi/full/10.1073/pnas.0600166103) simulates the SIRB with phage dynamics

- [Maity et al 2025](http://arxiv.org/abs/2503.21963) simulates the cholera driven by planktons

The key novelty in this current model includes:

- Combining both phages and planktons in understanding cholera outbreak

- Plankton-conferred protection on Vibrio against phages

vibriophage – minimal implementation + corrected README

This single document gives you:

Minimal R scripts to simulate five model variants (M0–M4) and generate synthetic observations

Data hooks to pull open cholera incidence and environmental (chlorophyll‑a) data

A corrected README with equations/notation consistent with your write‑up


## Model

The main model comprises three parts

- The basic SIR model
    - $S$: Susceptible human
    - $I$: Infected human
    - $R$: Recovered human
- Vibrio-phage interactions
    - $V_F$: free-living Vibrio cells
    - $V_Z$: zooplankton-associated Vibrio 
    - $P$: Phage
- Self-sufficient plankton dynamics:
    - $Z_F$: free-living Zooplanktons
    - $Z_V$: Vibrio-associated Zooplanktons
    - $X$: Phytoplanktons


$$
\begin{aligned}
\frac{dS}{dt} &= (b N + \omega R) - \delta S - \beta (V_F + V_Z) S  \\
\frac{dI}{dt} &= \beta (V_F + V_Z) S - (\gamma + \delta_I + \delta)I \\
\frac{dR}{dt} &= \gamma I - (\omega + \delta) R \\
\frac{dV_F}{dt} &= (\sigma I + \epsilon V_Z) - c \lambda V Z_F - \mu V_F P - \delta_V V \\
\frac{dV_Z}{dt} &= \rho Z_V - \mu g V_Z P - \epsilon Z_V - \delta_V V \\
\frac{dP}{dt} &= \tau \mu (V_F + g V_Z) P - \delta_P P\\
\frac{dZ_F}{dt} &= \eta \alpha X (Z_F + Z_V) - \lambda V Z_F - \delta_Z Z_F\\
\frac{dZ_V}{dt} &= \lambda V Z_F - \delta_Z Z_V\\
\frac{dX}{dt} &= r X (1 - \frac{X}{K}) - \alpha X (Z_V + Z_F)\\ 
\end{aligned}
$$

- $bN$ is the birth/immigration, where the total number of human population is $N = S+ I +R$
- $\omega R$ is the loss of immunity in recovered human
- $\delta S$ is the natural death/emmigration of human
- $\beta V S$ is the infection 


### Parameters

| Type | Parameter | Description  | Dimensions | Estimate  | Reference |
|:----:|:---------:|:------------------------------------------------------------|:-------------------------------:|:-----:|:-:|
SIR         | $b$       | Natural birth/immigration rate of human                  | day^-1^                         | $6.85 \times 10 ^{-5}$ | XX
SIR         | $\delta$    | Natural death/migration rate of human                    | day^-1^                         | $3.8 \times 10 ^{-5}$  | XX
SIR         | $\beta$     | Transmission rate via bacteria                           | day^-1^                         | $0.214$ | XX
SIR         | $\omega$    | Rate of immunity loss of recovered individuals           | day^-1^                         | $0.00092$ | XX
SIR         | $\gamma$    | Recovery rate of infected human                          | day^-1^                         | $0.2$  | XX
SIR         | $\delta_I$  | Induced death rate caused by infection                   | day^-1^                         | $0.013$ | XX
Vibriophage | $\sigma$    | Bacterial shedding rate by infected human                | cells L^-1^ day^-1^ per person  | $10$-$10^{4}$ | XX 
Vibriophage | $\rho$      | Bacterial shedding rate by infected zooplankton          | cells L^-1^ day^-1^ per mg      | XX | XX 
Vibriophage | $\mu$       | Encounter rate of Vibrio and phage                       | day^-1^                         | $1.4\times10^{-9}$ | XX
Vibriophage | $\delta_V$  | Death/washout rate of Vibrio in the reservoir            | day^-1^                         | $0.33$ | XX
Vibriophage | $\tau$      | Conversion coefficient from Vibrio to phage              | Virions per cell                | $80$ - $100$ | XX
Vibriophage | $\delta_P$  | Phage decay rate                                         | Virions per day                 | $0.5$-$7.9$ | XX
Plankton    | $\lambda$   | Encounter rate of Vibrio and zooplankton                 | day^-1^                         | $0.005$ - $0.1$ | XX
Plankton    | $c$         | Colonization efficient of bacteria in zooplankton        | cells (mg) L ^-1^               | $5 \times 10^{-7}$ | XX
Plankton    | $\eta$      | Conversion coefficient from phytoplankton to zooplankton | NA                              | $0.6$ | XX
Plankton    | $\alpha$    | Predation rate of zooplankton and phytoplankton          | day^-1^                         | $0.4$ | XX
Plankton    | $\delta_Z$  | Death rate of zooplankton                                | day^-1^                         | $0.06$ | XX
Plankton    | $r$         | Phytoplankton intrinsic growth rate                      | day^-1^                         | $0.5$ | XX
Plankton    | $K$         | Phytoplankton carrying capacity                          | mg L ^-1^                       | $1$ | XX


