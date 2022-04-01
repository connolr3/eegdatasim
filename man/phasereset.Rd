\name{phasereset}
\alias{phasereset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phase Reset
}
\description{
Generates EEG, assuming a Phase Reset Thoery Approach
}
\usage{
phasereset (frames, epochs, srate, minfr, maxfr , position, tjitter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frames}{number of signal frames per each trial}
  \item{epochs}{number of simulated trials}
  \item{srate}{sampling rate of simulated signal in Hz}
   \item{minfr}{minimum frequency of the sinusoid which is being reset}
  \item{maxfr}{maximum frequency of the sinusoid which is being reset}
    \item{postion}{position of the reset [in frames]; default: frames/2 => in the middle}
  \item{tjitter}{stdev of time jitter of the reset [in frames]; default: 0 => no jitter }
}
\details{
Function generates signal composed of a single sinusoid whose phase is
being reset. The frequency of the sinusoid is chosen randomly from a
specified range [minfr, maxfr]. The initial phase of the sinusoid is chosen randomly.
}
\value{
\item{signal }{simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials}
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
#myphasereset<-phasereset (200, 10, 250, 3, 9 , 115, 0)
#plot(myphasereset[1:200],type="l")

## The function is currently defined as
phasereset <-
  function (frames,
            epochs,
            srate,
            minfr,
            maxfr,
            position = frames / 2,
            tjitter = 0) {
    signal <- replicate(epochs * frames, 0)#generating empty wave
    for (trial in 1:epochs) {
      wavefr = runif(1, 0, 1) * (maxfr - minfr) + minfr
      initphase = runif(1, 0, 1) * 2 * pi
      pos = position + round(runif(1, 0, 1) * tjitter)
      for (i in 1:frames) {
        if (i < pos) {
          phase = i / srate * 2 * pi * wavefr + initphase
        }
        else{
          phase = (i - pos) / srate * 2 * pi * wavefr
        }
        signal[(trial - 1) * frames + i] = sin(phase)
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
}
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
