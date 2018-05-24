\name{Plotting for MMPPs}
\alias{MMPPplot}
\docType{package}
\title{Plotting function to visualise the results of the implemented Gibbs sampler}
\description{
This function provides visualisation of the results taken from either of the two implemented Gibbs samplers. The generated plot depicts the generated sample distribution of the latent Markov process \eqn{X_t} as well as the uncertainty of the results.

The blue line indicates the sample mode at each point in time while the thickness of the red lines indicates the amount of samples in each state. The uncertainty is calculated by \eqn{1-alpha_{t_i}} where \eqn{alpha_{t_i}} is the fraction of the sample mode of all samples at time $t_i$.

The plotting function relies on the R-packages 'ggplot2' and 'cowplot'.
}
\usage{
MMPPplot(Sampler_Output=NULL,
         title=" ",xaxis=" ",breaks=NULL,
         colour=NULL)       
}
\arguments{
  \item{Sampler_Output}{the output of one of the sampler functions "GibbsSampler" or "GibbsSampler_hierarchical"}
  \item{title}{an optional string stating the desired title of the plot}
  \item{xaxis}{an optional string stating the desired label for the x-axis}
  \item{breaks}{an optional dataframe with columns that replaces the x-axis ticks with customized tick labels. The first column indicates the position of the new ticks, while the second column contains strings with the new tick labels.}
  \item{colour}{an optional numeric vector that allows the colouring of specific regions in the data plot for better visual representation of specific time intervals. The vector indicates the position of the beginning/end of each colour interval.}
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
                                           B=1,N=5,messages=FALSE)

# Define appropriate new tick labels and colouring intervals
breaks <- data.frame(times=c(0,1800,
                          3600,5600,
                          7200,9000), 
                  names=c("14:00","14:30",
                          "15:00","15:30",
                          "16:00","16:30"))
colour <- c(0, 480, 1200, 2400, 
            2520, 3600, 4800, 6000, 
            7200, 7320, 8400, 9600)

example_plot <- MMPPplot(Sampler_Output=test_samples,
               title="Observations Imperial College Data",
               xaxis="time [hour]",
               breaks=breaks,
               colour=colour)
plot(example_plot)
}
