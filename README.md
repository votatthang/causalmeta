# CausalMeta

**Integration of aggregated data in causally interpretable meta-analysis by inverse weighting**

This repository contains the core R code and supporting files for a project that develops and evaluates methods for combining individual‑patient (IPD) and aggregated summary data in meta‑analysis under a causal inference framework. The approach relies on inverse‑weighting to produce effect estimates that respect causal assumptions while integrating multiple data levels.

---

## 📁 Repository contents

The repository consists primarily of a few scripts and accompanying text files:

- `DA_causal.R` — main analysis script implementing causal meta-analysis routines.
- `metaRE14.txt` — JAGS model specification files used in Bayesian analyses (loaded via `jags.model()` in R).

---

## 🛠️ Prerequisites

- **R** (version 4.0 or higher recommended).
- **JAGS** (Just Another Gibbs Sampler version 4.3.1) for Bayesian modeling (install from [mcmc-jags.sourceforge.io](https://sourceforge.net/projects/mcmc-jags/files/)).
- Required R packages (install from CRAN):

  ```r
  install.packages(c("nleqslv", "MASS", "rjags", "coda", "ggplot2"))
  ```

  These packages are loaded in the script:
  - `nleqslv`, `MASS`, `rjags`, `coda`, `ggplot2`.

---

## 🚀 Usage

1. **Clone the repository**

   ```sh
   git clone https://github.com/votatthang/causalmeta.git
   cd causalmeta
   ```

2. **Run analyses**

   Open `DA_causal.R` in R/RStudio and execute the script. 
   The script is self-contained and automatically loads the required data and functions from accompanying text files in the repository.

---

## 📝 Notes

- Data used by the script is not included; substitute your own datasets or adjust file paths accordingly.
- Results are printed directly to the R console. Inspect the script for details on the generated results.

---

## 📚 Citation

Please cite the associated manuscript or technical report when using these methods:

> *Integration of aggregated data in causally interpretable meta-analysis by inverse weighting*
> Vo, T. T., Le, T. T. K., Afach, S., & Vansteelandt, S.
> Submitted manuscript


---
