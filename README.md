# Replication Files for "Proximal Causal Inference for  Conditional Separable Effects"

This supplementary file contains replication codes for "Proximal Causal Inference for  Conditional Separable Effects" ([Park, Strensrud, Tchetgen Tchetgen, 2024](Link "PCSE")).

## Software and Packages

* Software: R version 4.2.1
* Packages: caret (version 6.0.93); earth (version 5.3.1); gam (version 1.22); gbm (version 2.1.8.1); glmnet (version 4.1.6); kernlab (version 0.9.32); MASS (version 7.3.58.2); matrix (version 1.4); mice (version 3.14.0); nnet (version 7.3.18); np (version 0.60.16); polspline (version 1.1.22); pracma (version 2.4.2); ranger (version 0.14.1); SuperLearner (version 2.0.28); xgboost (version 1.7.3.1)

## `PCSE_Simulation` Folder

* The `PCSE_Simulation` folder contains replication files for the simulation studies presented in Section 6 of the main paper and Section A9 of the supplementary material.

### Codes

* `0.MySL.R` contains functions for estimating nuisance functions using the Super Learner algorithm ([van der Laan, Polley, Hubbard, 2007](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html "SL")).

* `0.DGP.R` contains the data generating process.  

* `0.Functions\_PMMR.R` contains functions for estimating bridge functions based on the Proxy Maximum Moment Restriction (PMMR) method ([Mastouri et al., 2021](https://proceedings.mlr.press/v139/mastouri21a.html "PMMR")).  

* `1.Est\_PMMR\_Exp.R` and `1.Est\_PMMR\_Obs.R` 
	* These files replicate the simulation studies under experimental and observational settings, respectively.
	* Parallel computing is recommended for running these codes.
	* Each script returns an effect estimate using the proposed method.
	* Results are saved as `Result_PMMR_[aaa]_N[bbbb]_sdU[c]_B[ddddd].csv`, where:
	* `[aaa]` indicates whether the simulation data is generated from the observational or experimental setting.
	* `[bbbb]` indicates the number of observations in each simulated dataset.
	* `[c]` indicates whether an unmeasured confounder U is present or not.
	* `[ddddd]` indicates the random seed for generating a simulated dataset.

* `2.Est\_NoU\_Exp.R` and `2.Est\_NoU\_Obs.R`
	* These files replicate the simulation studies under experimental and observational settings, respectively.
	* Parallel computing is recommended for running these codes.
	* Each script returns an effect estimate under the no unmeasured confounding assumption.
	* Results are saved as `Result_Ign_[aaa]_N[bbbb]_sdU[c]_B[ddddd].csv`, where:
	* `[aaa]` indicates whether the simulation data is generated from the observational or experimental setting.
	* `[bbbb]` indicates the number of observations in each simulated dataset.
	* `[c]` indicates whether an unmeasured confounder U is present or not.
	* `[ddddd]` indicates the random seed for generating a simulated dataset.

* `3.Figure2_and_FigureA2.R` replicates Figure 2 in the main paper and Figure A2 in the supplementary material, visualizing effect estimates obtained from `1.Est\_PMMR\_Exp.R`, `1.Est\_PMMR\_Obs.R`, `2.Est\_NoU\_Exp.R`, and `2.Est\_NoU\_Obs.R`

* `4.Est\_PMMR\_NuisanceFt.R`
	* This file replicates the simulation study reported in the supplementary material.
	* Parallel computing is recommended.
	* This file computes the mean squared error of the PMMR nuisance function estimators.
	* Results are saved as `PSE_Nuisance_N[aa]_B[bbbbb].csv` files, where:
	* [aa] indicates the exponent of the number of observations in each simulated dataset.
	* [bbbbb] indicates the random seed for generating a simulated dataset.
	
* `5.FigureA3.R` replicates FigureA4 in the supplementary material, summarizing the convergence rates of the estimated nuisance functions from `4.Est\_PMMR\_NuisanceFt.R`.
	


## `PCSE_SimulatedData` Folder 

The `PCSE_SimulatedData` folder contains replication files that mimic the real-data analysis in Section 7 of the main paper. The Southwest Oncology Group (SWOG) data are not publicly available. Consequently, we provide replication code using simulated data that is compatible with the structure of the actual SWOG dataset. Because the simulated data may differ substantially from the true SWOG data, the resulting estimates may not match those reported in the main paper.

Interested readers can directly contact the study team of the SWOG trial (https://www.swog.org/) by requesting access to the dataset by requesting access via the provided email addresses at https://www.swog.org/about/contact-us. For a comprehensive understanding of the trial, please refer to  [Petrylak et al. (2004)](https://www.nejm.org/doi/full/10.1056/NEJMoa041318 "SWOG").

### Codes

* `0.MySL.R` contains functions for estimating nuisance functions using the Super Learner algorithm ([van der Laan, Polley, Hubbard, 2007](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html "SL")).

* `0.Functions\_PMMR.R` contains functions for estimating bridge functions based on the Proxy Maximum Moment Restriction (PMMR) method ([Mastouri et al., 2021](https://proceedings.mlr.press/v139/mastouri21a.html "PMMR")).  

* `1.Analysis.R`
	* This file replicates the data analysis using a simulated dataset.
	* Parallel computing is recommended for running this file.
	* Results are saved as `Result_PMMR_B[aaaaa].csv` and `Result_Ign_B[aaaaa].csv` files, where [aaaaa] indicates the random seed for generating a simulated dataset. These files correspond to the proposed approach and the no unmeasured confounding approach, respectively.

* `2.Table1.R` replicates Table 1 in the main paper, summarizing effect estimates from `1.Analysis.R`. As noted above, results may differ from those reported in the main paper due to the use of simulated data.

	

## References

Park, Stensrud, Tchetgen Tchetgen (2024). **Proximal Causal Inference for  Conditional Separable Effects**, _arXiv_ [[link](Link "PCSE")]

Mastouri et al. (2021). **Proximal causal learning with kernels: Two-stage estimation and moment restriction**, _International conference on machine learning_ [[link](https://proceedings.mlr.press/v139/mastouri21a.html "PMMR")]

van der Laan, Polley, Hubbard (2007). **Super learner**, _Statistical applications in genetics and molecular biology_ [[link](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html "SL")]

Petrylak et al. (2004). **Docetaxel and Estramustine Compared with Mitoxantrone and Prednisone for Advanced Refractory Prostate Cancer**, _The New England Journal of Medicine_ [[link](https://www.nejm.org/doi/full/10.1056/NEJMoa041318 "SWOG")]

