\name{est.sig.hat}
\alias{est.sig.hat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Hats
}
\description{
Estimates mean and standard deviation of noise part of signal
}
\usage{
est.sig.hats(data,peak_postion,buffer_pc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data on which to estimate from, should be averaged first}
  \item{peak_position}{Position of center of peak, if known, default is the index of max(data)}
  \item{buffer_pc}{Percentage of data on each side of peak_postion to exclude from calcualtions. Default is 0.3 which is very conservative as excludes 60pcof data total. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{1 }{Sample SD estimate of noise}
  \item{2 }{Sample mean estimate of noise}
%% ...
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
#mysignal<-noise(200,10,250)+peak(200,10,250,7,115)
#mysignal<-signal.averaging(mysignal,200,10)
#est.sig.hat(mysignal,115,0.15)

## The function is currently defined as


est.sig.hat<-function(data, peak_position=which.max(abs(data)),buffer_pc=0.3){
  lo<-peak_position-(buffer_pc*length(data))
  hi<-peak_position+(buffer_pc*length(data))
  buffer_range<-lo:hi
  reg_data<- data[-buffer_range]
  sig_hat<-sqrt(var(reg_data))
  normhat<-mean(reg_data)
  return(c(sig_hat,normhat))
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