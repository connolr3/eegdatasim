\name{signal_averaging1}
\alias{signal_averaging1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Signal Averaging (data, frames, epochs)
}
\description{
Averages Signal by transforming into matrix and obtaining column sums. Same result as signal_averaging2 function, just different methodology that may suit different size data.
}
\usage{
signal_averaging1(data,frames,trials)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector containing Signal To Average}
   \item{frames}{Frames in each trial}
    \item{epochs}{Number of Trials in vector}
}
\details{
Performs signal averging to drown out noise over a series of trials. By averaging a set of replicate measurements, the signal-to-noise ratio (SNR) will be increased,
}
\value{
\item{average_signal}{Averaged Signal of length (frames)}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rose Connolly connolr3@tcd.ie
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
https://en.wikipedia.org/wiki/Signal_averaging
}
\examples{
# mysignal<-noise(200,10,250)+7*peak(200,10,250,7,115)
#a_signal<-signal_averaging1(mysignal,200,10)
#plot_signal(a_signal,mark=115,title="Averaged Signal")


#reshape approach
signal_averaging1<-function(data,frames,epochs){
  m1 <- matrix(data, ncol=frames, byrow=TRUE)
  d1 <- as.data.frame(m1, stringsAsFactors=FALSE)
  summations<-colSums(d1)
  average_signal<-summations/epochs
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
