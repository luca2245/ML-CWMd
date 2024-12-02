# ML-CWMd

In this repository, you'll find an R implementation of ML-CWMd, a model rooted in the cluster-weighted models domain. It's designed for scenarios featuring:

- Hierarchical data
- Binary response
- Assumed presence of latent clusters among observations
- Multiple types of covariates: continuous, independent categorical, and dichotomous dependent covariates are all supported by the model.

## Methodology

- **U** $\rightarrow$ $p$-dimensional vector of continuous variables
- **V** $\rightarrow$ $q$-dimensional vector of categorical covariates
- **D** $\rightarrow$ $h$-dimensional vector of dichotomous covariates which may possess some degree of dependence
- **X** = (**U**,**V**,**D**) and Y (response variable) $\in$ $\mathbb{R}^{(p + q + h)}$ × {0,1} defined in a finite space $\boldsymbol{\Omega}$ which is assumed to be divided into $C$ clusters denoted as $\boldsymbol{\Omega}_1,\dots, \boldsymbol{\Omega}_C$
- Two-levels hierarchy, first-level observations $i$, for $i = 1, \dots, n_j$, are nested within groups $j$, for $j = 1, \dots, J$ 

- Joint probability across the clusters:
```math
p((\mathbf{x},\mathbf{y})| \boldsymbol{\theta}) = \sum_{c= 1}^{C}
p(\mathbf{y}|\mathbf{x},\boldsymbol{\xi}_{c})\phi(\mathbf{u}|
\boldsymbol{\mu}_{c},\boldsymbol{\Sigma}_c) \psi(\mathbf{v}|\boldsymbol{\lambda}_{c}) \zeta(\boldsymbol{d}| \boldsymbol{\Gamma}_c,\boldsymbol{\nu}_c)  w_{c} 
```
- $\boldsymbol{\theta}$ $\rightarrow$ vector containing all the parameters of the model 
- $w_c$ $\rightarrow$ indicates the proportion of observations within cluster c. 
- $\phi(\cdot|\boldsymbol{\mu}_c,\boldsymbol{\Sigma}_c)$ $\rightarrow$ multivariate normal density with cluster-wise different mean vectors $\boldsymbol{\mu}_c$ and covariance matrices $\boldsymbol{\Sigma}_c$
- $\psi(\cdot|\boldsymbol{\lambda}_c)$ $\rightarrow$ independent multinomial distributions with cluster-wise different parameter vectors $\boldsymbol{\lambda}_c$; 
- $\zeta(\cdot| \boldsymbol{\Gamma}_c,\boldsymbol{\nu}_c)$ $\rightarrow$ Ising model with cluster-wise different threshold vectors $\boldsymbol{\nu}_c$ and interaction matrices $\boldsymbol{\Gamma}_c$

A Shiny app is available at [https://mcdhmq-luca-caldera.shinyapps.io/app-ml-cwmd/](https://mcdhmq-luca-caldera.shinyapps.io/app-ml-cwmd/), allowing you to run the model on your data, visualize the results, and generate predictions.