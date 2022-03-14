\name{optimise_ERP}
\alias{optimise_ERP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
optimise_ERP
}
\description{
optimise uses OPTIM package to find parameters of ERP that minimise the SSE... It is equivalemt to fitting a model to a time series, however our model is based on trigonemtric functions. uses min_sse function as the function to minimise.
}
\usage{
optimise_ERP(dat,startingfrequency=0,starting_amp=1,mysr,pkcntr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{averaged signal from which to estimae from}
  \item{startingfrequency}{default=0, value for Optim to use as peak frequency starting value}
  \item{starting_amp}{default=1, value for Optim to use as peak amplitude starting value}
  \item{mysr}{Sampling Rate, in Hertz, used by experimenter}
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
# dat<-data.frame(x=erp_is,y=av_signal[erp_is])
#
# start_amp<-max(av_signal[erp_is])
# peak_center_estimate=which(av_signal==max(av_signal))
# pars<-optimise_ERP(dat,mysr=250,starting_amp = start_amp,pkcntr=peak_center_estimate)
#
#
# print(paste("estimating frequency as:", pars$par[1]))
# print(paste("Actual Frequency was:", 7))
# (paste("estimating amplitude as",pars$par[2]))
# print(paste("Actual Amplitude:", 3))


## The function is currently defined as
min_SSE<-function(data,par,srate,peakcenter){
  #  pred_values<-(1+cos(2*pi*par[1]*(data["x"]-par[2])))/2
  pred_values<-cos(((data["x"]-peakcenter)*2*pi*par[1])/srate)
  errs<-data["y"]-par[2]*pred_values
  sum((errs)^2)
}

#' @export
optimise_ERP<-function(dat,startingfrequency=0,starting_amp=1,mysr,pkcntr){
  result <- optim(par = c(startingfrequency, starting_amp), fn = min_SSE, data = dat,srate=mysr,peakcenter=pkcntr)
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
