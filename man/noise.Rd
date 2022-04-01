\name{noise}
\alias{noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Noise
}
\description{
Generates Classicaly simulated Noise from an EEG signal. Adapted from Matlab version.
Functioanlity is the same, however it avoids O(N^2) time complexity of Matlab version, as R is not suited well to for loops. Uses apply and matrices instead.
meanpower is read in at start in noise function, if not supplied by user.
}
\usage{
noise(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{frames}{number of signal frames per each trial}
  \item{epochs}{number of simulated trials}
  \item{srate}{sampling rate of simulated signal in Hz}
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
  \item{signal }{A vector of length frames*Epochs, representing (epochs) number of trials, each of length (frames), concatenated into one vector.}
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


fill_signals<-function(mysignal,frames, epochs, srate,meanpower){
  sumsig = 50#number of sinusoids from which each simulated signal is composed of
  signal <- matrix( rep(1:frames, sumsig), nrow=sumsig, byrow=TRUE )
  freq <- 4 * runif( sumsig, 0, 1)#generate random frequency for each sin
  freq <- cumsum(freq)#apply cumilitive function to make each sin have a higher frequency than the last
  freqamp <- meanpower[ pmin( ceiling(freq), rep(125,sumsig)) ] / meanpower[1]#generate ampltidue based on meanpower
  phase <- 2 * pi * runif(sumsig, 0 , 1)#generate random phase for each sin
  signal <- sin( signal * freq * ( 2 * pi / srate ) + phase ) * freqamp#create 50 sins, with consecutively higher freq for each sin 
  mysignal<-colSums( signal ) 
  
}

noise <- function(frames, epochs, srate, meanpower = NULL) {
  if (frames < 0) stop("frames cannot be less than 0")
  if (srate < 0) stop("srate cannot be less than 0")
  if (epochs < 0) stop("epochs cannot be less than 0")
  if( is.null(meanpower)) meanpower <- as.vector(R.matlab::readMat("meanpower.mat")$meanpower)
  signals<-matrix(0,frames,epochs)
  signals<-apply(signals,2,fill_signals,meanpower=meanpower,frames=frames,epochs=epochs,srate=srate)
  return(as.vector(signals))
}


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
