\name{Example data}
\alias{TestdataMMPP}
\alias{Testdataflows}
\docType{package}
\title{Example data to test the sampler}
\description{

These two datasets provide test examples to demonstrate the results of the sampler. They contain binned arrival data, the length of the binning intervals, the (assumed) number of process states, and appropriate prior hyperparameters. 

The first dataset, TestdataMMPP, contains artificial data generated from a Markov-modulated Poisson Process with two states, state decay rates \eqn{Q_ii}=0.01 for both states, and arrival rates \eqn{Lambda} with values 8 and 60 for good visible separation of the states.

The second dataset contains netflow arrival data collected from a personal computer by Clausen et. al (2017) <\url{https://github.com/hc2116/MMPPsampler/blob/master/Master_thesis_Henry.pdf}>, and can be used to test the implementation of the hierarchical MMPP model. The binning interval length is 5 seconds. It was concluded that the data is best described by a model with 5 states. For a more detailed description of the data, I point you to the given reference.

}
\author{
\packageAuthor{MMPPsampler}

Maintainer: \packageMaintainer{MMPPsampler}
}
\references{
Fearnhead, Paul, and Chris Sherlock. "An exact Gibbs sampler for the Markov-modulated Poisson process." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 68.5 (2006): 767-784 <\doi{10.1111/j.1467-9868.2006.00566.x}>

Clausen, Henry. "A Bayesian Approach to Human Behaviour Modelling in Computer Networks". Master's thesis, Imperial College London, <\url{https://github.com/hc2116/MMPPsampler/blob/master/Master_thesis_Henry.pdf}>
}
\seealso{
Source-code and more descriptions available under <\url{https://github.com/hc2116/MMPPsampler}>
}
\examples{
#Use the MMPP sample data included in the package to test the regular Gibbs sampler
data("TestdataMMPP")
Test <- TestdataMMPP
test_samples <- GibbsSampler(y_0T=Test$Bins,
                                    M=Test$M,
                                    Inter = Test$Inter,
                                    alpha_Gamma_rate = Test$alpha_Gamma_rate,
                                    alpha_Gamma_Q = Test$alpha_Gamma_Q,
                                    beta_Gamma_Q = Test$beta_Gamma_Q,
                                    beta_Gamma_rate = Test$beta_Gamma_rate,
                                    B=1,N=2,messages=FALSE)

test_plot <- MMPPplot(Sampler_Output = test_samples,
              title = "Example Plot")
plot(test_plot)

#Use the flow sample data included in the package to test the hierarchical model
data("Testdataflows")
Test <- Testdataflows
test_samples <- GibbsSampler_hierarchical(y_0T=Test$Bins,
                                           M=Test$M,
                                           Inter = Test$Inter,
                                           alpha_Gamma_rate = Test$alpha_Gamma_rate,
                                           alpha_Gamma_Q = Test$alpha_Gamma_Q,
                                           beta_Gamma_Q = Test$beta_Gamma_Q,
                                           beta_Gamma_rate = Test$beta_Gamma_rate,
                                           alpha_Gamma_Y=Test$alpha_Gamma_Z,
                                           beta_Gamma_Y=Test$beta_Gamma_Z,
                                           B=1,N=2,messages=FALSE)

}

