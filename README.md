# 📘 Joint Modeling of Sensitive Attribute and Observed Variable in Randomized Response Designs

This repository contains simulation code and results for the paper:

**“Joint Modeling of Sensitive Attribute and Observed Variable in Randomized Response Designs”**

Author: Shen-Ming Lee, Phuoc Loc Tran, Truong-Nhat Le, and Chin-Shang Li

## 📁 Files

| File Name                             | Description                                 |
|--------------------------------------|---------------------------------------------|
| `Supplement_Simulation_Study_01-02.R` | R script for simulation Case 1 & Case 2      |
| `Supplement_Simulation_Study_03.R`    | R script for simulation Case 3               |

## 🔬 Description

The simulations evaluate the performance of the proposed joint modeling method for estimating the relationship between a sensitive attribute (measured under RRT) and an observed variable. The results include comparisons of bias, standard deviation (SD), asymptotic standard error (ASE), and coverage probability (CP) across different conditions.

- **Case 1**: 
- **Case 2**: 
- **Case 3**: 

## ▶️ How to Run

1. Clone this repository or download the two `.R` files.
2. Open R or RStudio.
3. Install required packages (if not yet installed):

```r
install.packages(c("stats", "xtable"))
```

4. Run the script of interest:

```r
source("Supplement_Simulation_Study_01-02.R")
# or
source("Supplement_Simulation_Study_03.R")
```

## 📊 Output

Each script will:

- Run 1000 simulations (default setting)
- Estimate model parameters under different conditions
- Output tables summarizing:
  - Parameter bias
  - Empirical SD
  - Estimated ASE
  - Coverage probability (CP)
  - LRT test
