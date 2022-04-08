\name{optimise_ERP}
\alias{optimise_ERP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
optimise_ERP
}
\description{
optimise uses OPTIM package to find parameters of ERP that minimise the SSE... Usually, using OPTIM in R is equivalemt to fitting a model to a time series. However our model is based on trigonemtric functions.

optimise_ERP uses min_sse function as the function to minimise.
It uses BFGS method to optimise. The gradient of the SSE is supplied also, which speeds up computation (gradient is the derivative of the SSE function).

Altough optime_ERP estimate both amplitude and frequency..... the optimisation is only one-dimensional as the amplitude can be formed in terms of the freqeuncy. 
}
\usage{
optimise_ERP(dat,startingfrequency=0,starting_amp=1,mysr,pkcntr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{x values from data we are fitting}
  \item{y}{y values from data we are fitting}
  \item{mysr}{Sampling Rate, in Hertz, used and knwon by experimenter}
  \item{pkcntr}{Index of peak center, if unknown please estimate with the which(max(dat))}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{par }{The best set of parameters found. par[1] = freq, par[2]=amplitude}
\item{value }{The value of fn corresponding to par.}
\item{counts }{A two-element integer vector giving the number of calls to fn and gr respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.}
\item{convergence }{The value of fn corresponding to par. An integer code. 0 indicates successful completion (which is always the case for "SANN" and "Brent"). Possible error codes are

1
indicates that the iteration limit maxit had been reached.

10
indicates degeneracy of the Nelder--Mead simplex.

51
indicates a warning from the "L-BFGS-B" method; see component message for further details.

52
indicates an error from the "L-BFGS-B" method; see component message for further details.}
\item{message }{Any additional info returned by optimiser, may be Null}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rose Connolly connolr3@tcd.ie
}
\note{
Please see: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/optim
for R documentation on Optim Function
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# mysignal<-noise(200,10,250)+3*peak(200,10,250,7,115)
# av_signal<-signal_averaging1(mysignal,200,10)
# hats<-est_sig_hat(av_signal)
# stan_signal<-(av_signal-hats[2])/hats[1]
# erp_is<-find_ERP_range(stan_signal,cutoff=1.6)
#
#
# peak_center_estimate=which(av_signal==max(av_signal))
# pars<-optimise_ERP(erp_is,av_signal[erp_is],mysr=250,peak_center_estimate)
#
#
# print(paste("estimating frequency as:", pars$par[1]))
# print(paste("Actual Frequency was:", 7))
# (paste("estimating amplitude as",pars$par[2]))
# print(paste("Actual Amplitude:", 3))


## The function is currently defined as
gr_min_SSE <- function(par,x,y,srate,peakcenter){
  u <- ((x-peakcenter)*2*pi)/srate
  z <- cos( u * par[1] )
  alphaest <- sum( z * y) / sum(z*z)
  ans<-4*pi*alphaest*(x-peakcenter)*(y-alphaest*z)*(sin(u * par[1]))
  return(sum(ans/srate))
}

min_SSE<-function(par,x,y,srate,peakcenter){
  z <- cos(((x-peakcenter)*2*pi*par[1])/srate)
  alphaest <- sum(z * y) / sum(z*z)
  #pred_values<-cos(((data["x"]-peakcenter)*2*pi*par[1])/srate)
  errs<- y - alphaest * z #pred_values
  sum(errs^2)
}

#' @export
optimise_ERP<-function(x,y,mysr,pkcntr){
  result <- optim(par = 1, fn = min_SSE, gr=gr_min_SSE, x=x, y=y,srate=mysr,peakcenter=pkcntr, method="BFGS")
  z <- cos(((x-pkcntr)*2*pi*result$par[1])/mysr )
  alphaest <- sum(z * y) / sum(z*z)
  result$par <- c((result$par),alphaest)
  return(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
}