#  Reduction of Autocorrelation Times in Lattice Path Integral Quantum Monte Carlo via Direct Sampling of the Truncated Exponential Distribution

Emanuel Casiano-Diaz, Kipton Barros, Ying Wai Li, Adrian Del Maestro

[arXiv:XXXX.XXXXX]

### Abstract
 In Monte Carlo simulations, proposed configurations are accepted or rejected according to an acceptance ratio, which depends on an underlying probability distribution and an *a priori* sampling probability. By carefully selecting the probability distribution from which random variates are sampled,  simulations can be made more efficient, by virtue of an autocorrelation time reduction. In this paper, we illustrate how to directly sample random variates from a two dimensional truncated exponential distribution. We show that our direct sampling approach converges vaster to the target distribution as compared to rejection sampling. The direct sampling of one and two dimensional truncated exponential distributions is then applied to a recent Path Integral Monte Carlo (PIMC) algorithm for the simulation of Bose-Hubbard lattice models at zero temperature. The new sampling method results in improved acceptance ratios and reduced autocorrelation times of estimators without any wall clock run time increase, providing an effective speed up of the simulation.


### Description
This repository includes links, code, scripts, and data to generate the plots in the above paper.

### Requirements

The figures in this project were generated from the data files in the [processed_data](https://github.com/DelMaestroGroup/papers-code-truncExponSampling/tree/master/processed_data) directory of this repository with the various `.py` scripts found in the [src](https://github.com/DelMaestroGroup/papers-code-truncExponSampling/tree/master/src) directory. Figures were generated using the `.ipynb` notebook files contained there.

The raw data for Figures 6,7,8 was generated via a new path integral Monte Carlo algorithm for the ground state of bosonic lattice models: [pigsfli](https://github.com/DelMaestroGroup/pigsfli).

### Support
The creation of these materials was supported in part by the National Science Foundation under Award Nos. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991&HistoricalAwards=false) and [DMR-2041995](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2041995&HistoricalAwards=false).

[<img width="100px" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

This work was supported by the Laboratory Directed Research and Development Early Career Research Funding of Los Alamos National Laboratory (LANL) under project number 20210662ECR. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration of U.S. Department of Energy (Contract No.
89233218CNA000001). A.~D. was supported by the U.S.  Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Number DE-SC0022311.

### Figures

#### Figure 1: One dimensional truncated exponential distribution
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/simple_truncexpon_histogram_benchmark.svg" width="400px">

#### Figure 2: Kolmogorov-Smirnov (KS) test 
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/0.100000_1.500000_-2.000000_simpleTruncexpon_ksTest.svg" width="400px">

#### Figure 3: Two dimensional truncated exponential distribution
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/joint_truncexpon_histogram_benchmark.svg" width="400px">

#### Figure 4: Cumulative averages vs number of samples
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/rvs_convergence.svg" width="400px">

#### Figure 6: System size dependence of integrated autocorrelation time
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/pigsfli_tau_vs_L_critical.svg" width="400px">

#### Figure 7: Projection length dependence of integrated autocorrelation time
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/pigsfli_tau_vs_beta_critical.svg" width="400px">

#### Figure 8: Comparison of PIMC wall clock times with rejection and direct sampling
<img src="https://github.com/DelMaestroGroup/papers-code-truncExponSampling/blob/master/figures/wall_times.svg" width="400px">

Figures are relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

