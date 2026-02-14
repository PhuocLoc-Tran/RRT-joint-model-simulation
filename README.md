# 📘 Conditional distributional framework for joint modeling of sensitive attribute and an observed Variable under randomized response designs

This repository contains the simulation and empirical analysis scripts accompanying the paper:

**“Conditional distributional framework for joint modeling of sensitive attribute and an observed Variable under randomized response designs”**
Shen-Ming Lee, Phuoc Loc Tran, Truong-Nhat Le, and Chin-Shang Li
(JRSSC revision)

The scripts reproduce the simulation studies and empirical analyses reported in the manuscript.

---

# 📂 Repository Structure

| File Name                                   | Description                                                                |
| ------------------------------------------- | -------------------------------------------------------------------------- |
| `Supplement_Simulation_Study_01-02.R`       | Simulation code for Case 1 and Case 2 (Section 4 of the paper)             |
| `Supplement_Simulation_Study_03.R`          | Simulation code for Case 3 (Section 4 of the paper)                        |
| `Empirical_Analysis_TSCS.R`                 | Analysis code for the 2012 Taiwan Social Change Survey example (Section 5) |

---

# 🔬 Purpose of the Simulation Study

The simulation experiments evaluate the finite-sample performance of the proposed joint likelihood framework for modeling:

$$ P(Y, Z \mid \boldsymbol{X}) = P(Y \mid \boldsymbol{X}) P(Z \mid Y, \boldsymbol{X}).$$

where:

* $Y$: latent sensitive binary variable collected under the unrelated-question RRT,
* $Z$: observed binary response variable,
* $\boldsymbol{X}$: covariates.

The simulations assess:

* Bias of parameter estimates
* Empirical standard deviation (SD)
* Average asymptotic standard error (ASE)
* Coverage probability (CP) of Wald confidence intervals
* Likelihood ratio test (LRT) performance for testing $H_0: \alpha_0 = \alpha_1$.

All simulation designs correspond exactly to Section 4 of the manuscript.

---

# 🧪 Simulation Design

Each script implements Monte Carlo experiments under the following settings:

* Sample size: $n \in \{1000, 2000\}$
* RRT design parameters: $p \in \{0.5, 0.7\},\ c \in \{0.25, 0.5\}$
* Number of replications: 1000 (default).

## Case 1

* $X_1$: continuous
* $X_2$: ordinal

## Case 2

* $X_1$: binary
* $X_2$: ordinal

## Case 3

* $X_1$: binary
* $X_2$: binary

The true parameter values used in each scenario match those reported in the manuscript tables.

---

# ⚙️ Structure of the R Scripts

Each simulation script follows the same workflow:

1. **Generate covariates**
   Simulate $X_1$, $X_2$ according to the specified case.

2. **Generate latent sensitive variable $Y$**
   Using logistic model:
   $$P(Y=1|\boldsymbol{X}) = H(\boldsymbol{\beta}^\top \boldsymbol{X}).$$

3. **Generate observed response $Z$**
   Using conditional model:
   $$P(Z=1|Y=y,\boldsymbol{X}) = H(\boldsymbol{\alpha}_y^\top \boldsymbol{X}).$$

4. **Apply unrelated-question RRT mechanism**
   Generate randomized response $Y^*$.

5. **Estimate parameters using EM algorithm**

   * E-step: compute posterior expectations of $Y$.
   * M-step: update logistic regression parameters.

6. **Compute performance metrics**

   * Bias
   * SD
   * ASE
   * CP
   * LRT rejection frequency.

Extensive inline comments are provided within each script explaining these steps.

---

# ▶️ How to Run the Code

## Step 1: Clone or Download

Download the `.R` scripts to your local machine.

## Step 2: Open R or RStudio

## Step 3: Install Required Packages

## Step 4: Run a Simulation Script


---

# 📊 Output

Each script:

* Runs 1000 Monte Carlo replications
* Stores parameter estimates across replications
* Produces summary tables including:

  * Bias
  * Empirical SD
  * Average ASE
  * Coverage probability (95% CI)
  * LRT rejection rate

The reported values correspond to the simulation tables and figures in the manuscript.

---

# 📈 Empirical Example (Section 5)

The empirical analysis is based on the 2012 Taiwan Social Change Survey (TSCS).

Due to data governance policies, the raw data cannot be redistributed in this repository.
However, researchers may apply for access through the Survey Research Data Archive (SRDA), Academia Sinica.

The analysis scripts included here:

* Define the binary variables $Y$, $Z$
* Specify covariates $X_1$, $X_2$.
* Fit the joint model
* Compute parameter estimates and LRT statistics.

Researchers who obtain access to the TSCS data through SRDA can directly reproduce all empirical results.

---

# 🔁 Reproducibility Notes

* All simulation settings match those described in Section 4 of the paper.
* Random seeds are set to ensure replicability.
* Results may vary slightly due to Monte Carlo variability.
* The EM algorithm implementation follows the likelihood formulation described in Section 3.
