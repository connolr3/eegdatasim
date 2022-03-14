\name{peak}
\alias{peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak
}
\description{
Generates Classicaly simulated Peak (ERP) from an EEG signal
}
\usage{
peak(frames,epochs,srate,peakfr,position = frames / 2,tjitter = 0,wave = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frames}{number of signal frames per each trial}
  \item{epochs}{number of simulated trials}
  \item{srate}{sampling rate of simulated signal in Hz}
  \item{peakfr}{frequency of sinusoid whos half of the cycle is taken to form the peak}
  \item{position}{position of the peak [in frames]; default: frames/2 => in the middle}
  \item{tjitter}{stdev of time jitter of the peak; default: 0 => no jitter}
  \item{wave}{if defined the signal is composed not from a peak, but complete sinusoid. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Signal }{Vector of length frames*epochs containing peak part of a simulated EEG signal}

}
\references{
https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
}
\author{
Adapted from Yeung N,Bogacz R, Holroyd C, Nieuwenhuis S, Cohen J
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#mypeak<-peak(200,10,250,7,115)
#plot(mypeak[1:200])


## The function is currently defined as

peak <-
  function(frames,
           epochs,
           srate,
           peakfr,
           position = frames / 2,
           tjitter = 0,
           wave = NULL) {
    mypeak <- c(1, 5, 4, 9, 0)
    mypeak
    signal <- replicate(epochs * frames, 0)
    for (trial in 1:epochs) {
      pos = position + round(runif(1, 0, 1) * tjitter)

      for (i in 1:frames) {
        phase = (i - pos) / srate * 2 * pi * peakfr
        if ((is.null(wave) == FALSE) ||
            (phase < pi / 2 &&
             phase > -pi / 2)) {
          #if wave | (phase < pi/2 & phase > -pi/2)
          signal[(trial - 1) * frames + i] = cos(phase)
        }
      }
    }
    return(signal)
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
