\name{Makinen}
\alias{Makinen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makinen
}
\description{
Generates Signal as per Makinen et. al 2005 paper
}
\usage{
Makinen (frames, epochs, srate, position)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{frames}{number of signal frames per each trial}
  \item{epochs}{number of simulated trials}
  \item{srate}{sampling rate of simulated signal in Hz}

   \item{position=frames/2}{position of the reset [in frames]; default: in the middle}
}
\details{
Generate simulated data by summing 4 phase reset sinusoids with freqencies chosen randomly from range 4-16Hz. Has the same parameters as phasereset except the parameters describing the frequency range.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{signal }{simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
#mymak = Makinen (200, 1, 250, 115)
#plot (mymak)

## The function is currently defined as
Makinen <- function(frames, epochs, srate, position = frames / 2) {
  signal <- replicate(epochs * frames, 0)#generating empty wave
  #print(length(signal))
  for (i in 1:4) {
    #repeat for 4 sinusoids
    new_signal <- phasereset(frames, epochs, srate, 4, 16, position)
    #print(length(new_signal))
    signal = signal + new_signal#add new sin to summation
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
}