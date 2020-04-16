# PartIII-Astrostatistics-2020
Webpage for Astrostatistics Course, Lent Term 2020  
Part III Mathematics / Part III Astrophysics

**Example Sheet 4 has been uploaded!**  

Lent Term
Monday, Wednesday & Friday at 12 noon. CMS **Meeting Room 5**.

First Class: Friday 17 Jan 2020

Office Hours: Tuesdays @ 1:30pm or by appointment  
Statistical Laboratory  
CMS Pavilion D, Office 1.07  

Recommended Texts:  
(F&B and Ivezić are freely available as downloadable electronic PDFs through the Cambridge University Library iDiscover website.)

* F&B = Feigelson & Babu. "Modern Statistical Methods for Astronomy"  
* Ivezić = Ivezić, Connolly, VanderPlas & Gray. "Statistics, Data Mining, and Machine Learning in Astronomy"  
* MacKay = David MacKay. "Information Theory, Inference, and Learning Algorithms"  
(Free text: http://www.inference.org.uk/itila/book.html)  
* Givens & Hoeting. "Computational Statistics" (Free through Cambridge Library iDiscover)  
* Bishop, "Pattern Recognition and Machine Learning"  
(Free text: https://www.microsoft.com/en-us/research/people/cmbishop/#!prml-book)  

**Lecture 01 - 17 Jan 2020**  
Introduction to the Course

**Lecture 02 - 20 Jan 2020**  
Introduction to Astronomical Data & Case Studies  

**Lecture 03 - 22 Jan 2020**  
Probability foundations - F&B Ch 2, Ivezić Ch 3  
properties of random variables, conditional probability, Bayes' Theorem  
multivariate continuous rvs, transformation of rvs  

**Lecture 04 - 24 Jan 2020**  
Probability foundations -  
multivariate pdfs, covariance, propagation of error, limit theorems  
discussion of estimators  

**Lecture 05 - 27 Jan 2020**  
statistical inference - F&B Ch 3, Ivezić Ch 4  
Frequentist properties of estimators: unbiasedness, consistency, efficiency,  
mean square error, bias-variance tradeoff  

**Lecture 06 - 29 Jan 2020**  
deriving estimators: method of moments, least squares  
maximum likelihood, properties of MLE, Fisher information, Cramer-Rao Bound  

**Week 3**  
**Lecture 07 - 31 Jan 2020** - slides posted  
* Multi-parameter MLE, Fisher matrix
* Local Distance Ladder  
* Estimating the absolute magnitude distribution of supernovae with measurement error  

**Lecture 08 - 03 Feb 2020** - slides posted  
* Calibrating the absolute magnitude distribution of supernovae with heteroskedastic measurement error 
* Fitting a power law (Pareto distribution) with observational selection effects
* Quantifying uncertainty

**Lecture 09 - 05 Feb 2020**  
* interpetation of frequentist confidence intervals
* Quantifying uncertainty with bootstrap
* Example: estimating a Pareto distribution with selection effects  
* begin regression  

**Week 4** - Statistical Modelling & Bayesian Inference in Astrophysics  
**Lecture 10 - 07 Feb 2020**
* Review standard regression methods: OLS, WLS, GLS
* Probabilistic forward (generative) modelling with latent variables
* Regression with heteroskedastic (x,y)-errors and intrinsic dispersion
* Kelly et al. "Some Aspects of Measurement Error in Linear Regression of Astronomical Data"  
* Comparison of regression methods
* Statistical Modelling wisdom

**Lecture 11 - 10 Feb 2020**  
* Bayesian Modelling and Inference
* Frequentist confidence intervals vs Bayesian credible intervals
* Simple Gaussian example w/ conjugate prior  
  - Limit of flat prior
  - Large Sample limit  

**Lecture 12 - 12 Feb 2020**  
* Credible intervals: central vs. highest posterior density (HPD)
* Example: Bailer Jones et al. "Estimating Distances from Parallax"
* Use of "noninformative" vs physically-motivated priors  
* Frequency evaluation of Bayesian estimators  

**Example Class 1 - 13 Feb 2020**  
* Local Hubble Constant estimation with Supernova and Cepheid data  

**Week 5** - Statistical Computation  
**Lecture 13 - 14 Feb 2020**  
* Intuitive discussion of local Hubble constant estimation  
* Finish parallax discussion
* More Bayesian Inference: Gaussian with unknown mean and variance  
* Sufficient statistics, Likelihood principle  
* Conjugate Priors, Large Sample Limits, Asymptotic Normality  
* Monte Carlo sampling, Direct simulation  

**Lecture 14 - 17 Feb 2020**
* Continue Monte Carlo Sampling  
* Gaussian with unknown mean and variance, noninformative priors
* Direct sampling
* Importance Sampling and Self-Normalised Importance Weights  
* Application to Bayesian parallax problem
* Application: Patel et al. 2017, 2018 - Weighing the Milky Way Galaxy  

**Lecture 15 - 19 Feb 2020**  
* Application: Patel et al. 2017, 2018 - Weighing the Milky Way Galaxy  
* Markov Chain Monte Carlo: 1D Metropolis Algorithm  

**Week 6** - More Markov Chain Monte Carlo  
**Lecture 16 - 21 Feb 2020**  
* 1D Metropolis applied to Parallax example  
* d-dimensional Metropolis Algorithm  
* Multivariate Gaussian draws using Cholesky factors  
* Gaussian with unknown mean and variance example  

**Lecture 17 - 24 Feb 2020**
* Supernova cosmology example with 4D Metropolis  

**Lecture 18 - 26 Feb 2020**  
* continue SN cosmology 4D Metropolis  
* Assessing convergence with multiple chains and Gelman-Rubin  
* Metropolis-Hastings algorithms  
* Gibbs sampling: correlated 2D Gaussian example  

**Week 7**  
**Lecture 19 - 28 Feb 2020**  
* Metropolis-within-Gibbs  
* Supernova cosmology example: dark energy equation-of-state  
* Assessing mixing: autocorrelation time, thinning  
* Efficiency comparison of different MCMC algorithms  

**Lecture 20 - 02 Mar 2020**  
* Begin sketch of MCMC theory  
* demo of convergence of MCMC to stationary dist'n  
* aperiodic, irreducible, ergodic chains  
* stationary distributions  

**Lecture 21 - 04 Mar 2020**
* Finish sketch of MCMC theory  
* detailed balance, reversibility  
* Metropolis(-Hastings) algorithm respects DB  
* Mixed Metropolis/Gibbs samplers and parameter blocking  
* MCMC in practice  
* Begin Gaussian Processes  

**Week 8**  
**Lecture 22 - 06 Mar 2020**
* Gaussian Processes and Astronomical Time Series Analysis  
* Gravitationally Lensed SN example  

**Lecture 23 - 09 Mar 2020**  
* Hierarchical Bayesian models in astrophysics  
* Probabilistic graphical models and Conditional Independence  
* "Normal-Normal" multi-level model  

**Lecture 24 - 11 Mar 2020**  
* Gibbs sampling for hierarchical models  
* Shrinkage and partial pooling  
* Bayesian model comparison, Bayes factors, and evidence  
