# Bayesian R-Pareto Mixture for Spatial Extremes

This repository contains the code for my thesis on modelling extreme rainfall using a **Bayesian mixture of \(\mathcal{R}\)-Pareto processes**. The goal is to jointly model different physical drivers of extreme precipitation (frontal vs. convective events) in a single probabilistic framework.

## Key Ideas

- Models multiple extremal regimes **jointly**, not separately.
- Learns:
  - spatial range parameters (Brown–Resnick semivariogram),
  - latent membership (frontal vs. convective),
  - mixture weights.
- Uses an \(L_1\)-based \(\mathcal{R}\)-Pareto likelihood; extendable to other \(L_p\).
- Validated through simulation and applied to real rainfall from Tampa Bay, Florida.

## Repository Structure

```text
├─ code                 
│  ├─ application_final.Rmd # code for application study
│  └─ simulation_final.Rmd  # code for simulation
└─ README.md

@misc{rpareto-mixture-sevt,
  author  = {Haitao Gao},
  title   = {Bayesian Mixture of \mathcal{R}-Pareto Processes for Spatial Rainfall Extremes},
  year    = {2025},
  note    = {Thesis codebase},
  url     = {https://github.com/accuracy-maker/gradute_thesis.git}
}
