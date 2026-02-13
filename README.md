# KAMA: a knockoff-augmented meta-analysis framework for genome-wide association studies

**KAMA** (**K**nockoff-**A**ugmented **M**eta-**A**nalysis) is a meta-analysis framework for simultaneously identifying risk variants and pinpointing the specific studies in which they are significant. It supports various data types, including individual-level data and summary statistics, by integrating knockoff statistics derived from each study. KAMA performs meta-analysis with rigorous False Discovery Rate (FDR) control, enhancing the power to identify risk variants. Moreover, it localizes the identified signals to specific studies, providing insights into the sources of genetic risk and enhancing the interpretability of meta-analysis findings (e.g., in multi-ancestry and multi-design GWAS).

This package provides functions to calculate knockoff statistics, perform meta-analysis with FDR control, and localize signals.

![workflow](/figure/KAMA_workflow.png)

---

## Repository Structure

- `R/KAMA.R`: Defines the main meta-analysis and localization functions.
    - `KAMA_meta_analysis`: Computes $q$-values for FDR control.
    - `KAMA_localization`: Localizes signals to specific studies.
- `R/stats.R`: Contains statistical utility functions.
    - `MK.statistic`: Computes $\kappa$ and $\tau$ statistics from importance scores.
- `R/data.R`: Documentation for example datasets.
- `data/study.assoc_stats.rda`: Example association statistics ($p$-values) from 4 studies.
- `data/pops.sum_res.rda`: Example summary statistics (Z-scores) and LD matrices for 4 populations, formatted for GhostKnockoff generation.
- `figure/`: Stores the workflow diagram.
- `man/`: Documentation files (`.Rd`) for all functions.

---

## Installation

You can install the development version of KAMA from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install KAMA from GitHub
devtools::install_github("JianYin-Gloria/KAMA")
```

**Required Dependencies**:
- R (>= 4.0.0)
- `dplyr`, `magrittr`, `stats`, `tibble`, `tidyr`

---

## Application 1: Meta-Analysis with Marginal P-values
This application illustrates the core workflow of KAMA using marginal association statistics. KAMA can directly use marginal $p$-values for both original variants and their knockoff replicates, regardless of whether they were derived from individual-level genotypes or summary statistics.

### **Step 0**: Load Data
The file (included in the package) `data("study.assoc_stats")` contains simulated association statistics ($p$-values) regarding genomic regions across four independent studies: `study1`, `study2`, `study3`, `study4`.
It serves as input for calculating knockoff statistics ($\kappa$ and $\tau$) and performing subsequent meta-analysis and localization.

```r
library(KAMA)
library(dplyr)

data("study.assoc_stats")
```

The data is a named list of 4 data frames. Each data frame corresponds to a study and contains:
- `chr`: Integer. The chromosome number.
- `start`, `end`: Integer. The start and end position of the genomic region.
- `actual_start`, `actual_end`: Integer. The position of the first and last variant in the region.
- `p_value`: Numeric. The $p$-value from the original association test.
- `p_value_kc1` ... `p_value_kc5`: Numeric. $p$-values from 5 knockoff replicates.

### Step 1: Compute Multiple Knockoff Statistics
First, transform $p$-values into importance scores ($T$) and calculate multiple knockoff statistics ($\kappa$ and $\tau$).

```r
M <- 5 # Number of knockoffs

# Compute Kappa and Tau for each study
kappa_list <- list()
tau_list <- list()

for (study in names(study.assoc_stats)) {
  df <- study.assoc_stats[[study]]
  
  # Transform p-values to importance scores (larger is better)
  T_0 <- -log10(df$p_value)
  T_k <- -log10(as.matrix(df[, paste0("p_value_kc", 1:M)]))
  
  # Compute Statistics
  stats <- MK.statistic(T_0, T_k)
  
  kappa_list[[study]] <- stats[, "kappa"]
  tau_list[[study]]   <- stats[, "tau"]
}

kappa_set <- do.call(cbind, kappa_list)
tau_set   <- do.call(cbind, tau_list)

# Optional: Get genetic info for annotation
genetic_info <- study.assoc_stats[[1]] %>% select(chr, start, end, actual_start, actual_end)
```

#### Arguments
- **`T_0`**: A numeric vector or a one-column matrix representing the importance scores of the **original** features (variables).
- **`T_k`**: A numeric matrix representing the importance scores of the **knockoff** features. It must have the same number of rows as `T_0` and `M` columns (where `M` is the number of knockoff replicates).
- **`method`**: A character string specifying the aggregation method for $\tau$. Currently supports `'median'`. Default is `'median'`.

#### Outputs
A matrix with two columns:
- **`kappa`**: The index of the original (denoted as 0) or knockoff feature (1 to `M`) that has the largest importance score. 
- **`tau`**: The difference between the largest importance score and the median of the remaining importance scores.

### Step 2: Run Meta-Analysis
Apply the KAMA framework to identify significant variants across all studies.

```r
meta_results <- KAMA_meta_analysis(
  kappa_set = kappa_set, 
  tau_set = tau_set, 
  M = M, 
  genetic_info = genetic_info
)

# Filter for significant results (e.g., q < 0.1)
sig_results <- meta_results %>% filter(qvalue <= 0.1)
head(sig_results)
```

#### Arguments
- **`kappa_set`**: A numeric matrix of $\kappa$ statistics. Rows correspond to variables, and columns correspond to studies. **Note:** The order of variables (rows) and studies (columns) must strictly match that of `tau_set`.
- **`tau_set`**: A numeric matrix of $\tau$ statistics. Rows correspond to variables, and columns correspond to studies. **Note:** The order of variables (rows) and studies (columns) must strictly match that of `kappa_set`.
- **`M`**: Integer. The number of knockoffs used. Must be >= max(kappa_set).
- **`genetic_info`**: Optional data frame or matrix containing additional information (e.g., SNP IDs, positions) to be appended to the results.
- **`seed`**: Integer. Random seed for reproducibility (default 111).

#### Outputs
A data frame containing:
- **`genetic_info`** (if provided): The appended metadata.
- **`kappa_set`** columns: The input $\kappa$ statistics.
- **`tau_set`** columns: The input $\tau$ statistics.
- **`qvalue`**: KAMA's $q$-value for meta-analysis FDR control.

### Step 3: Localization
Localize the identified signals to specific studies.

```r
loc_results <- KAMA_localization(
  kappa_set = kappa_set,
  tau_set = tau_set,
  M = M,
  genetic_info = genetic_info,
  q = 0.1,
  study_label = names(study.assoc_stats)
)

print(loc_results)
```

#### Arguments
- **`kappa_set`**, **`tau_set`**, **`M`**, **`genetic_info`**, **`seed`**: Same as `KAMA_meta_analysis`.
- **`q`**: Numeric. The FDR control level (default 0.1).
- **`study_label`**: Character vector. Custom labels for the studies (columns). Defaults to "study1", "study2", etc.

#### Outputs
A data frame of the **significant** variables containing:
- **`genetic_info`** (if provided): The appended metadata.
- **`kappa_set`**, **`tau_set`**: The knockoff statistics for the significant variables.
- **`qvalue`**: KAMA's $q$-value for meta-analysis FDR control.
- **`localization`**: A character column specifying the studies/cohorts where the variable is significant.

---

## Application 2: Meta-Analysis from Summary Statistics & LD
This workflow assumes you have standard GWAS summary statistics (Z-scores) and LD matrices. It uses the `GhostKnockoff` package (external) to generate knockoff scores without individual-level data.

### **Step 0**: Load Data
The file `data("pops.sum_res")` contains summary statistics and LD matrices for 233 variants across four independent populations (`AFR`, `EAS`, `EUR`, `SAS`).

```r
library(GhostKnockoff)
data("pops.sum_res")
```
The data is a named list of 4 elements (AFR, EAS, EUR, SAS). Each element contains:
- `sum_stats`: A data frame with 233 rows and 9 columns, including:
  - `chr`, `pos`, `SNP`: Variant identifiers.
  - `A1`, `A2`: Reference and alternative alleles.
  - `Beta`: Effect size estimate.
  - `SE`: Standard error.
  - `Zscore`: Z-statistic (Beta / SE).
- `cor_mat`: A 233 x 233 symmetric numeric matrix representing the LD correlation structure.

### Step 1: Generate Knockoff Statistics
Use the provided `pops.sum_res` dataset which contains Z-scores and LD matrices.

```r
# Note: This example requires the 'GhostKnockoff' package to generate statistics.
# library(GhostKnockoff)
# data(pops.sum_res)

# Parameters
M <- 5 # number of knockoffs
n.study <- 2000 # sample size per population
kappa_list <- list()
tau_list <- list()

# 2. Process each population
for (pop in names(pops.sum_res)) {
  sum_stats <- pops.sum_res[[pop]]$sum_stats
  cor_mat   <- pops.sum_res[[pop]]$cor_mat

  # Generate knockoffs (Set seed for reproducibility)
  set.seed(99)
  fit.prelim <- GhostKnockoff.prelim(cor_mat, M = M, method = 'asdp', max.size = 500)

  # Fit knockoff statistics
  GK.stat <- GhostKnockoff.fit(
    Zscore_0 = as.matrix(sum_stats$Zscore),
    n.study = n.study,
    fit.prelim = fit.prelim,
    gamma = 1
  )

  # Calculate KAMA statistics (Kappa and Tau)
  stats <- MK.statistic(GK.stat$T_0, GK.stat$T_k)
  kappa_list[[pop]] <- stats[, "kappa"]
  tau_list[[pop]]   <- stats[, "tau"]
}

# 3. Prepare input matrices
kappa_set <- do.call(cbind, kappa_list)
tau_set   <- do.call(cbind, tau_list)

# Add prefixes to distinguish columns
colnames(kappa_set) <- paste0("kappa_", names(pops.sum_res))
colnames(tau_set)   <- paste0("tau_", names(pops.sum_res))

# Extract genetic information from the first study
genetic_info <- pops.sum_res[[1]]$sum_stats[, c("chr", "SNP", "pos", "A1", "A2")]
```

#### Outputs
- **`kappa_set`**: A matrix combining $\kappa$ statistics from all populations.
- **`tau_set`**: A matrix combining $\tau$ statistics from all populations.

### Step 2: Run Meta-analysis and Localization
Once `kappa_set` and `tau_set` are obtained, apply the standard KAMA pipeline.

```r
# Run Meta-analysis
meta_res <- KAMA_meta_analysis(kappa_set, tau_set, M = M, genetic_info = genetic_info)

# Run Localization
loc_res  <- KAMA_localization(kappa_set, tau_set, M = M, q = 0.1,
                              genetic_info = genetic_info,
                              study_label = names(pops.sum_res))

print(loc_res)
```

#### Outputs
See the "Outputs" section in **Application 1** for details on the return values of `KAMA_meta_analysis` and `KAMA_localization`.

---

## Contact

For questions or issues, please open an GitHub issue or contact the maintainer:
**Jian Yin** (jian.yin@my.cityu.edu.hk).
