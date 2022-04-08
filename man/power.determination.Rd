\name{power.determination}
\alias{power.determination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Determination 
}
\description{
Obtains an estimate of bernoulli power based on monte carlo estimation.

Obtains this power in terms of ERP frequency estimation, and ERP ampltide estimation.
}
\usage{
power.determination(0.1,7,5,200,250,40,50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{accuracy_window}{percentage used to caluclate range acceptable estimate is found in}
  \item{freq}{expected peak frequency }
    \item{amp}{expected peak amplitude}
      \item{frames}{Frames in each trial}
        \item{srate}{Maximum number of trials you wish to calculate p for, default = 40}
          \item{srate}{number of trials for wish to avaergae over to estimate p for each trial, default = 100}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{results}{dataframe containing ps for ERP frequency estimations in first column, and ps for ERP amplitude estimations in second.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rose Connolly connolr3@tcd.ie
}
\note{
Has a significant run-time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#NOT RUN
# my.ps<-power.determination(0.1,7,5,200,250,30,50)
#plot(my.ps$my_frequency_ps,main="p estimates for 1 to 30 trials averaged over 50 simulations")

## The function is currently defined as

power.determination <-function(accuracy_window,freq,amp,frames,srate,maxtrial = 40,averagingN=100)
{
  set.seed("123")
  my_frequency_ps <- numeric(length=maxtrial)
  my_amplitude_ps <- numeric(length=maxtrial)
  meanpower <- R.matlab::readMat("meanpower.mat")$meanpower
  for (N in 1:maxtrial)
  {
    freq_ones <- 0
    amp_ones <- 0
    for (j in 1:averagingN)
    {
      #create sample data
      mysignal <-noise(frames, N, 250, meanpower=meanpower) + amp * peak(frames, N, srate, freq)
      
      #prepare signal: average and standardise
      my_averaged_signal <- signal_averaging(mysignal, frames, N)
      hats <- est_sig_hat(my_averaged_signal, frames / 2)
      standata <- (my_averaged_signal - hats[2]) / hats[1]
      
      #estimate peak range
      mypeak_range <- find_ERP_range(standata, 1.7)
      
      #prepare data and starting values for optim
      yis <- my_averaged_signal[mypeak_range]
      peak_center_estimate = which(my_averaged_signal == max(my_averaged_signal))
      
      pars <-optimise_ERP(mypeak_range, yis ,mysr = srate,pkcntr = peak_center_estimate)
      if (abs(freq - pars$par[1]) <= freq * accuracy_window) 
        freq_ones=freq_ones+ 1
      if (abs(amp - pars$par[2]) <= amp * accuracy_window) 
        amp_ones=amp_ones+ 1
    }
    
    cat("\n Completed (no. trials): ",N)
    freq_p <- freq_ones / averagingN
    amp_p <- amp_ones / averagingN
    my_frequency_ps[N]<-freq_p
    my_amplitude_ps[N]<-amp_p
  }
  results<-data.frame(my_frequency_ps,my_amplitude_ps)
  return((results))
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
