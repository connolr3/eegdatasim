\name{find_ERP_range}
\alias{find_ERP_range}
roxygen2::roxygenise()
\title{
find_ERP_range
}
\description{
Estimates ERP range in a signal
}
\usage{
find_ERP_range(data,cutoff)
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
#mysignal<-noise(200,10,250)+peak(200,10,250,7,115)
#av_signal<-signal_averaging1(mysignal,200,10)
#hats<-est_sig_hat(av_signal)
#stan_signal<-(av_signal-hats[2])/hats[1]
#erp_is<-find_ERP_range(stan_signal,cutoff=1.6)
#mydata<-data.frame(av_signal,col="black")
#col<-rep("black",200)
#col[erp_is]="red"
#plot(av_signal,col=col)


find_ERP_range<-function(data,cutoff=2){
  index<-which.max( abs(data) )
  i<-index+1
  while(abs(data[i][1])>cutoff & i<length(data))  #goright
  {
    index<-c(index,i)
    i<-i+1
  }
  i<-index[1]-1
  while(abs(data[i][1])>cutoff & i<length(data))#goleft
  {

    index<-c(i,index)
    i<-i-1
  }
  return(index)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~eeg }
\keyword{ ~erp }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
