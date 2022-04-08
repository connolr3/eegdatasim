\name{plot.signal}
\alias{plot.signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Signal
}
\description{
Quickfire way to Plot EEG noise, peak, or signal
}
\usage{
plot.signal(data,mark,title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector to Plot}
  \item{mark}{Optional int, add horixontal line marker}
  \item{title}{Optional string, add title}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%
%%  \item{comp2 }{Description of 'comp2'}
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
# mysignal<-noise(200,10,250)+7*peak(200,10,250,7,115)
#firstsignal<-mysignal[1:200]
#plot.signal(firstsignal,mark=115,title="Classicaly generated Signal with Peak")

plot.signal <- function(data, mark = NULL, title = "EEG Signal") {
  plot(data, main = title, type = "l")
  # this if statement and corresponding code is applicable for phase resetting theory only
  if (is.null(mark) == FALSE) {
    abline(v = mark, col = "red")
    text(mark + 5, 0.25, paste(mark), srt = 90, col = "red")
  }
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