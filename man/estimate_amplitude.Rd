\name{estimate_amplitude}
\alias{estimate_amplitude}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimate_amplitude
}
\description{
Estimates Amplitude of Signal as Max value
}
\usage{
estimate_amplitude(averaged_signal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{averaged_signal}{Averaged Signal for function to be performed on
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{ }{Maximum value of data serving as estimate for peak amplitude}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#mysignal<-noise(200,10,250)+peak(200,10,250,7,115)
#mysignal<-signal_averaging1(mysignal,200,10)
#estimate_amplitude(mysignal)

## The function is currently defined as
estimate_amplitude<-function(averaged_signal){
  return(max(averaged_signal))
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
