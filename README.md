Bayesian Modeling with Spatial Curvature Processes
================

<## Aritra Halder, Sudipto Banerjee, Dipak K. Dey>

This GitHub repository contains necessary code/scripts required to reproduce the results in, "Bayesian Modeling with Spatial Curvature Processes". The required sub-routines to reproduce the analysis can be found in https://github.com/arh926/spWombling. The repository also outlines a workflow using a simulated dataset.

## Data

The manuscript uses three datasets, (i) Boston Housing Data (ii) Meuse River Data (iii) Temperatures in the Northeastern US which are available in .Rda format in the arh926/spWombling/data folder. They can also be independently accessed through R-packages `spData` (https://cran.r-project.org/web/packages/spData/index.html) and `spBayes` (https://cran.r-project.org/web/packages/spBayes/index.html).


## Code
The code contained in this repository is instrumental in producing the tables and plots in the manuscript. We provide brief descriptions:

`sim_sp_final.R`: This is the script that, upon replication, is responsible for producing the Tables S1 and S2 in the online Supplement that documents simulations that assess the accuracy of estimated gradients and curvature. Script for generating plots for gradients, curvature and other differential geometric constructs are also included in the script. Although only Pattern 1 is highlighted, Pattern 2 is commented out and can be run if required. The script concludes by computing wombling measures for curves shown in Pattern 1 producing Table 1 of manuscript.

`meuse-spatial.R`: This produces the application results for Meuse River data (second part of Section 6). This contains the full application, including computation of wombling measures.

`boston-spatial.R`: This produces the application results for Boston Housing Data (first part of Section 6). This concludes with outlining curves of interest  for the dataset which can then be used as an input in `bayes_cwomb.R` (in the package) or `cwomb-riemann.R` (provided) to compute wombling measures shown in the paper.

`netemp-spatial.R`: This produces the application results for Temperatures in the Northeastern US Data (Supplement). This concludes with outlining curves of interest  for the dataset which can then be used as an input in bayes_cwomb.R or cwomb-riemann.R to compute wombling measures shown in the paper.

`sim-jasa.R`: Contains additional simulation investigating effect of the spatial field's variance on the width of highest posterior density intervals for gradienst and curvature.

The subroutines in https://github.com/arh926/spWombling/ R-package contain detailed descriptions. The README file also shows workflow for estimating gradients and wombling measures under Pattern 1.

## Instructions for use
For running each of the above scripts successfully, installing the R-package from https://github.com/arh926/spWombling/ is advised. Instructions for installation are provided in the README.


