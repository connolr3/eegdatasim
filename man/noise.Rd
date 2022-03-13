\name{noise}
\alias{noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Noise
}
\description{
Generates Classicaly simulated Noise from an EEG signal
}
\usage{
noise(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frames }{Frames in one Trial}
  \item{Epochs }{Trials conducted in experiment}
  \item{srate }{Sampling rate in Hz}
}
\details{
The value of the first parameters describing the number of samples is computed by multiplying the duration of the noise by the sampling frequency, i.e. if we wanted 0.8 seconds of noise, with a srate of 250Hz.... choose frames = 200
0.8 * 250 = 200
The function generates a vector containing the samples.

Noise is generated such that its power spectrum matches the power spectrum of human EEG
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{signal }{A vector of length frames*Epochs, representing different trials of equal length each }
%% ...
}
\references{
https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
}
\author{
Rose Connolly connolr3@tcd.ie
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


#noise <- function(frames, epochs, srate) {
 # sumsig = 50#number of sinusoids from which each simulated signal is composed of
  ##signal <- matrix(0,1,epochs*frames)
  #signal <- replicate(epochs * frames, 0)
  #for (trial in 1:epochs) {
   # freq = 0
  #  signal_range <- c(((trial - 1) * frames + 1):(trial * frames))
   # for (i in 1:sumsig) {
    #  freq = freq + (4 * runif(1, 0, 1))
     # freqamp = meanpower[min (ceiling(freq), 125)] / meanpower[1]#meanpower created up above from read-in file
      phase <- runif(1, 0, 1) * 2 * pi
  #    signal[signal_range] <-
   #     signal[signal_range] + sin(c(1:frames) / srate * 2 * pi * freq + phase) * freqamp
#
#    }
#  }
#  signal
#  return(signal)
#}

#mynoise<-noise(200,10,250)
#plot(mynoise[1:200])


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
