\name{find.ERP.range}
\alias{find.ERP.range}
\title{
find.ERP.range
}
\description{
Estimates ERP range in a signal
}
\usage{
find.ERP.range(data,cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data }{Vector (should be an averaged, standardised signal) in which to find ERP}
  \item{cutoff=2 }{Minimum value of point to be included in ERP.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{index}{Vector of indexes containing estimated ERP}
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
find.ERP.range<-function(data,cutoff=2){
  z <- abs(data)
  z <- z - cutoff
  index<-which.max( z )
  
  zi <- z > 0
  
  left.side <- rev( zi[1:index] )
  t <- which( left.side == FALSE )
  low <- index - min(t) + 1
  
  right.side <- zi[index:length(data)]
  t <- which( right.side == FALSE )
  high <- index + min(t) - 1

  return( low:high )
}
#mysignal<-noise(200,10,250)+peak(200,10,250,7,115)
#av_signal<-signal.averaging(mysignal,200,10)
#hats<-est.sig.hat(av_signal)
#stan_signal<-(av_signal-hats[2])/hats[1]
#erp_is<-find.ERP.range(stan_signal,cutoff=1.6)
#mydata<-data.frame(av_signal,col="black")
#col<-rep("black",200)
#col[erp_is]="red"
#plot(av_signal,col=col)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
