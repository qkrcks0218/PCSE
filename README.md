# Replication Files for "Proximal Causal Inference for  Conditional Separable Effects"

This supplementary file contains replication codes for "Proximal Causal Inference for  Conditional Separable Effects".

## Software and Packages

* Software: R version 4.2.1
* Packages: caret (version 6.0.93); earth (version 5.3.1); gam (version 1.22); gbm (version 2.1.8.1); glmnet (version 4.1.6); kernlab (version 0.9.32); MASS (version 7.3.58.2); matrix (version 1.4); mice (version 3.14.0); nnet (version 7.3.18); np (version 0.60.16); polspline (version 1.1.22); pracma (version 2.4.2); ranger (version 0.14.1); SuperLearner (version 2.0.28); xgboost (version 1.7.3.1)

## Common Code Files

* 0.MySL.R contains functions used for estimating nuisance functions based on the Superlearner algorithm ([van der Laan, Polley, Hubbard, 2007](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html "SL")).
* 0.Functions\_PMMR.R contains functions used for estimating bridge functions based on the Proxy Maximum Moment Restriction (PMMR) method ([Mastouri et al., 2021](https://proceedings.mlr.press/v139/mastouri21a.html "PMMR")).  

### Code 

* 1.Est\_PMMR.R
	* This file replicates the simulation study. Parallel computing is recommended. 
	* The results are saved as "Result_PMMR_[aaa]_N[bbbb]_B[ccccc].csv" files in "Result" folder.
	* [aaa] indicates whether the simulation data is generated from the observational or experimental setting.
	* [bbbb] indicates the number of observations in each simulated dataset.
	* [ccccc] indicates the random seed for generating a simulated dataset.
* 2.EffectCalculate.R
	* This file calculates the true target effects. Parallel computing is recommended.
	* The results are saved as "Effect_B[aaaaa].csv" files in "Effect" folder.
	* [aaaaa] indicates the random seed for generating a simulated dataset. 
* 3.Summary\_Merge.R
	* This file aggregates csv files in "Result" and "Effect" folder, respectively, and generates a PNG file summarizing the simulation result. 

## References

Mastouri et al. (2021). **Proximal causal learning with kernels: Two-stage estimation and moment restriction**, _International conference on machine learning_ [[link](https://proceedings.mlr.press/v139/mastouri21a.html "PMMR")]

van der Laan, Polley, Hubbard (2007). **Super learner**, _Statistical applications in genetics and molecular biology_ [[link](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html "SL")]
