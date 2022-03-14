\name{signal_averaging2}
\alias{signal_averaging2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Signal Averaging 2
}
\description{
Averages Signal by iterating over entire vector. More efficient for smaller data. Same result as signal_averaging1 function, just different methodology that may suit different size data.
}
\usage{
signal_averaging2(data,frames,epochs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector containing Signal To Average}
   \item{frames}{Frames in each trial}
    \item{epochs}{Number of Trials in vector}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# mysignal<-noise(200,10,250)+7*peak(200,10,250,7,115)
#a_signal<-signal_averaging2(mysignal,200,10)
#plot_signal(a_signal,mark=115,title="Averaged Signal")

## The function is currently defined as
#approach by summation
signal_averaging2<-function(data,frames,epochs){
  average_signal<-rep(0,frames)
  range<- c(1:frames)
  for (i in c(1:epochs)) {
    average_signal<-average_signal+data[range]
    range<-range+frames
  }
  return(average_signal/epochs)
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
