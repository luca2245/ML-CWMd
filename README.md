# ML-CWMd

In this repository, you'll find an R implementation of ML-CWMd, a novel methodology that integrates and extends cluster-weighted models (CWM) with logistic mixed-effects models (GLMM). The ML-CWMd is able to identify latent subpopulations, capture hospital-specific effects, and manage the conditional dependencies resulting from dichotomous covariates, incorporating the Ising model within the CWM framework. Additionally, it properly accounts for the dependencies introduced by the nested data structure. 

It is specially designed for complex data scenarios involving:

- **Hierarchical structures**
- **Binary responses**
- **Latent clusters** within observations
- **Diverse covariates**: continuous, independent categorical, and dependent dichotomous variables

## Model Overview

### Key Components
The model works with:

- **U**: A $p$-dimensional vector of continuous variables
- **V**: A $q$-dimensional vector of categorical covariates
- **D**: An $h$-dimensional vector of dichotomous covariates, which may exhibit dependence

Together, these covariates form $\mathbf{X} = (\mathbf{U}, \mathbf{V}, \mathbf{D})$, and the response variable $Y \in \mathbb{R}^{(p + q + h)} \times \{0,1\}$. The data space, $\boldsymbol{\Omega}$, is assumed to consist of $C$ latent clusters: $\boldsymbol{\Omega}_1, \dots, \boldsymbol{\Omega}_C$.

### Hierarchical Structure
The model incorporates a two-level hierarchy:

1. **First Level**: Observations $i$, for $i = 1, \dots, n_j$
2. **Second Level**: Groups $j$, for $j = 1, \dots, J$

### Joint Probability Across Clusters
The model estimates the joint probability as follows:

```math
p((\mathbf{x}, \mathbf{y})| \boldsymbol{\theta}) = \sum_{c=1}^{C} p(\mathbf{y}|\mathbf{x}, \boldsymbol{\xi}_{c}) \phi(\mathbf{u}| \boldsymbol{\mu}_{c}, \boldsymbol{\Sigma}_c) \psi(\mathbf{v}|\boldsymbol{\lambda}_{c}) \zeta(\mathbf{d}| \boldsymbol{\Gamma}_c, \boldsymbol{\nu}_c) w_{c}
```

- **$\boldsymbol{\theta}$**: Vector of all model parameters
- **$w_c$**: Proportion of observations in cluster $c$
- **$\phi(\cdot| \boldsymbol{\mu}_c, \boldsymbol{\Sigma}_c)$**: Multivariate normal density with cluster-specific means ($\boldsymbol{\mu}_c$) and covariances ($\boldsymbol{\Sigma}_c$)
- **$\psi(\cdot| \boldsymbol{\lambda}_c)$**: Independent multinomial distributions with cluster-specific parameters ($\boldsymbol{\lambda}_c$)
- **$\zeta(\cdot| \boldsymbol{\Gamma}_c, \boldsymbol{\nu}_c)$**: Ising model with cluster-specific thresholds ($\boldsymbol{\nu}_c$) and interactions ($\boldsymbol{\Gamma}_c$)

---

## Getting Started

### Run the Model Locally
You can test the model directly in R using the script **`run_algorithm.R`**. Example datasets are provided in the "Example Datasets" folder, featuring the two distinct data-generating processes (DGPs) considered in the paper:

- **`data1.csv`**: Three latent clusters
- **`data2.csv`**: Two latent clusters

#### Part 1: Quick Testing with Provided Datasets
1. Choose a dataset and initialize the model with:
   - **Random initialization**, or
   - **K-means initialization**
2. Specify the number of latent clusters.
3. Evaluate the prediction performance using the corresponding test datasets:
   - **`test_data1.csv`** for `data1.csv`
   - **`test_data2.csv`** for `data2.csv`

#### Part 2: Full Cluster Selection Procedure
Follow the complete procedure outlined in the paper to identify the optimal number of latent clusters:

1. **Define a Range for `C`**: Specify the range of possible cluster values.
2. **Run the Algorithm for Each `C`**:
   - For each value of `C`, run the algorithm 30 times with random initialization and select the iteration with the highest log-likelihood.
   - If all random initializations fail, switch to k-means initialization.
3. **Compute BIC**: Calculate the Bayesian Information Criterion (BIC) for the best iteration of each `C`.
4. **Choose the Optimal `C`**: Select the number of clusters corresponding to the lowest BIC value.

### Explore the Model with a Shiny App
A **Shiny app** is available to test the ML-CWMd model, visualize results, and generate predictions. Access the app here: [https://mcdhmq-luca-caldera.shinyapps.io/app-ml-cwmd/](https://mcdhmq-luca-caldera.shinyapps.io/app-ml-cwmd/). You can use the example datasets provided in the "Example Datasets" folder to test the model within the app.

---