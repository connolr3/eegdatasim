
#See: https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator


usethis::use_package("R.matlab")
usethis::use_package("optimr")


#' @export
plot_signal <- function(data, mark = NULL, title = "EEG Signal") {
  plot(data, main = title, type = "l")
  # this if statement and corresponding code is applicable for phase resetting theory only
  if (is.null(mark) == FALSE) {
    abline(v = mark, col = "red")
    text(mark + 5, 0.25, paste(mark), srt = 90, col = "red")
  }
}

#' @export
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


noise_y <- function(frames, epochs, srate, meanpower = NULL) {
  if (frames < 0) stop("frames cannot be less than 0")
  if (srate < 0) stop("srate cannot be less than 0")
  if (epochs < 0) stop("epochs cannot be less than 0")
  if( is.null(meanpower)) meanpower <- as.vector(R.matlab::readMat("meanpower.mat")$meanpower)
  sumsig = 50#number of sinusoids from which each simulated signal is composed of
  #signal <- matrix(0,1,epochs*frames)
  signal <- replicate(epochs * frames, 0)
  signaln<-array(rep(0, sumsig*epochs*frames), c(sumsig, epochs, frames))
  for (trial in 1:epochs) {
    #freq = 0
    signal_range <- ((trial - 1) * frames + 1):(trial * frames)
    S <- matrix( rep(1:frames, sumsig), nrow=sumsig, byrow=TRUE )
    f <- 4 * runif( sumsig, 0, 1)
    f <- cumsum(f)
    fa <- meanpower[ pmin( ceiling(f), rep(125,sumsig)) ] / meanpower[1]
    phi <- 2 * pi * runif(sumsig, 0 , 1)
    S <- sin( S * f * ( 2 * pi / srate ) + phi ) * fa
    signal[ signal_range ] <- colSums( S )    
  }
  signal
  return(signal)
}






#' @export
peak <- function(frames, epochs,srate,peakfr,position = frames / 2,tjitter = 0,wave = NULL) {
  if (frames < 0) stop("frames cannot be less than 0")
  if (srate < 0) stop("srate cannot be less than 0")
  if (epochs < 0) stop("epochs cannot be less than 0")
    mypeak <- c(1, 5, 4, 9, 0)
    mypeak
    signal <- replicate(epochs * frames, 0)
    for (trial in 1:epochs) {
      pos = position + round(runif(1, 0, 1) * tjitter)

      for (i in 1:frames) {
        phase = (i - pos) / srate * 2 * pi * peakfr
        if ((is.null(wave) == FALSE) ||
            (phase < pi / 2 &&
             phase > -pi / 2)) {
          #if wave | (phase < pi/2 & phase > -pi/2)
          signal[(trial - 1) * frames + i] = cos(phase)
        }
      }
    }
    return(signal)
  }


#' @export
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


#' @export
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


#' @export
Makinen1a <- function() {
  trials = 30
  mysignal = Makinen (400, trials, 1000, 175)#30 trials
  my_new_signal <-
    matrix(mysignal, nrow = 400, ncol = 30)#reshape into matrix
  my_df <-
    as.data.frame(t(my_new_signal))#convert -> df to use a matlab plot

  par(mfrow = c(3, 1))
  #plot1
  matplot(t(my_df), type = "l", ylab = "EEG")
  #plot2
  signalmean <- lapply(my_df[1:400], FUN = mean)
  plot(
    c(1:400),
    signalmean,
    type = "l",
    col = "blue",
    xlab = "",
    ylab = "ERP"
  )
  #plot3
  variance <- lapply(my_df[1:400], FUN = var)
  plot(c(1:400),
       variance,
       type = "l",
       col = "blue",
       xlab = "")
  par(mfrow = c(1, 1))#reset
}






#approach by summation

#' @export
signal_averaging_x<-function(data,frames,epochs){
  average_signal<-rep(0,frames)
  range <- c(1:frames)
  for (i in 1:epochs) {
    average_signal<-average_signal+data[range]
    range<-range+frames
  }
  return(average_signal/epochs)
}

signal_averaging<-function(data,frames,epochs){
  a <- matrix( data, nrow=frames, ncol=epochs )
  return(rowMeans(a))
}

#' @export
estimate_amplitude<-function(averaged_signal){
  return(max(averaged_signal))
}

#' @export
est_sig_hat_print<-function(data, peak_position=which.max(data),buffer_pc=0.2){
  #data - *averaged* signal we are working with
  #peak_position - center of peak
  #buffer_pc - % of either side of peak we want to exclude from sig hat estimation .... the higher the more conservative

  lo<-peak_position-(buffer_pc*length(data))
  hi<-peak_position+(buffer_pc*length(data))
  buffer_range<-lo:hi
  reg_data<- data[-buffer_range]
  sig_hat<-sd(reg_data)
  normhat<-mean(reg_data)

  #plotting and printing part.... no actual functionality
  par(mfrow=c(2,1))
  plot(data,main=paste("Identified Noise part of data ",buffer_pc*100,"% either side of peak"))
  abline(v=c(lo,hi),col="red")
  hist(reg_data,main="Histogram of Noise part of Signal")
  print(paste("SD of data is: ",sqrt(var(reg_data))))
  print(paste("Mean of data is: ",mean(reg_data)))
  print(paste("Data is of length",length(data)," and we are calculating noise sd & mean based on values outside the range",lo," to ",hi,". This is within a range of ",buffer_pc*100,"% of the peak center at ",peak_position))
  par(mfrow=c(1,1))

  return(c(sig_hat,normhat))
}


#' @export
est_sig_hat<-function(data, peak_position=which.max(abs(data)),buffer_pc=0.3){
  lo<-peak_position-(buffer_pc*length(data))
  hi<-peak_position+(buffer_pc*length(data))
  buffer_range<-floor(lo):floor(hi)
  reg_data<- data[-buffer_range]
  sig_hat<-sd(reg_data)
  normhat<-mean(reg_data)
  return(c(sig_hat,normhat))
}

#' @export
find_ERP_range<-function(data,cutoff=2){
  z <- abs(data)
  z <- z - cutoff
  index<-which.max( z )
  
  zi <- z > 0
  
  left_side <- rev( zi[1:index] )
  t <- which( left_side == FALSE )
  low <- index - min(t) + 1
  
  right_side <- zi[index:length(data)]
  t <- which( right_side == FALSE )
  high <- index + min(t) - 1

  return( low:high )
}

#plots signal with erp highlighted in red
plot_erp<-function(signal,erp_range){
  x_values<-c(1:length(signal))
  point_colour<-replicate(length(signal),"black")
  point_colour[erp_range]<-"red"
  mydf<-data.frame(x_values,signal,point_colour)
  plot(x_values,signal,col=point_colour,main="Signal with ERP identified")
}

gr_min_SSE <- function(par,x,y,srate,peakcenter){
  u <- ((x-peakcenter)*2*pi)/srate
  z <- cos( u * par[1] )
  alphaest <- sum( z * y) / sum(z*z)
  r <- x - alphaest * z
  v <- alphaest * -sin( u * par[1] ) * u
  return( 2 * sum( r * v ) )
}

gr_min_SSE2 <- function(par,x,y,srate,peakcenter){
  u <- ((x-peakcenter)*2*pi)/srate
  z <- cos( u * par[1] )
  alphaest <- sum( z * y) / sum(z*z)
  ans<-4*pi*alphaest*(x-peakcenter)*(y-alphaest*z)*(sin(u * par[1]))
  return(sum(ans/srate))
}

min_SSE<-function(par,x,y,srate,peakcenter){
  z <- cos(((x-peakcenter)*2*pi*par[1])/srate)
  alphaest <- sum(z * y) / sum(z*z)
  #pred_values<-cos(((data["x"]-peakcenter)*2*pi*par[1])/srate)
  errs<- y - alphaest * z #pred_values
  sum(errs^2)
}

#' @export
optimise_ERP<-function(x,y,mysr,pkcntr){
  result <- optim(par = 1, fn = min_SSE, gr=gr_min_SSE2, x=x, y=y,srate=mysr,peakcenter=pkcntr, method="BFGS")
  z <- cos(((x-pkcntr)*2*pi*result$par[1])/mysr )
  alphaest <- sum(z * y) / sum(z*z)
  result$par <- c((result$par),alphaest)
  return(result)
}


power_determination <-function(accuracy_window,freq,amp,frames,srate,freq_power = TRUE,
                               amp_power = FALSE,maxtrial = 100,averagingN=4)
{
  meanpower<- R.matlab::readMat("meanpower.mat")$meanpower
  meanpower <- as.vector(meanpower)
  set.seed("123")
  my_ps <- numeric()
  for (N in 1:maxtrial)
  {
    ones <- 0
    for (j in 1:averagingN)
    {
      #create sample data
      mysignal <-noise(frames, N, 250,meanpower) + amp * peak(frames, N, srate, freq)
      
      #prepare signal: average and standardise
      my_averaged_signal <- signal_averaging(mysignal, frames, N)
      hats <- est_sig_hat(my_averaged_signal, frames / 2)
      standata <- (my_averaged_signal - hats[2]) / hats[1]
      
      #estimate peak range
      mypeak_range <- find_ERP_range(standata, 1.7)
      
      #prepare data and starting values for optim
      yis <- my_averaged_signal[mypeak_range]
      peak_center_estimate = which.max(my_averaged_signal)
      
      pars <-optimise_ERP(mypeak_range, yis ,mysr = srate,pkcntr = peak_center_estimate)
      if (abs(freq - pars$par[1]) <= freq * accuracy_window) {
        ones=ones+ 1
      }
      
    }
    cat("\n Completed (no. trials): ",N)
    p <- ones / averagingN
    my_ps[N]<-p
  }
  return(my_ps)
}


power_determination_amp <-function(accuracy_window,freq,amp,frames,srate,maxtrial = 100,averagingN=4)
{
  set.seed("123")
  my_ps <- c()
  meanpower <- R.matlab::readMat("meanpower.mat")$meanpower
  for (N in 1:maxtrial)
  {
    ones <- 0
    for (j in 1:averagingN)
    {
      #create sample data
      mysignal <-noise(frames, N, 250, meanpower=meanpower) + amp * peak(frames, N, srate, freq)
      
      #prepare signal: average and standardise
      my_averaged_signal <- signal_averaging(mysignal, frames, N)
      hats <- est_sig_hat(my_averaged_signal, frames / 2)
      standata <- (my_averaged_signal - hats[2]) / hats[1]
      
      #estimate peak range
      mypeak_range <- find_ERP_range(standata, 1.7)
      
      #prepare data and starting values for optim
      yis <- my_averaged_signal[mypeak_range]
      peak_center_estimate = which(my_averaged_signal == max(my_averaged_signal))
      
      pars <-optimise_ERP(mypeak_range, yis ,mysr = srate,pkcntr = peak_center_estimate)
      if (abs(amp - pars$par[2]) <= amp * accuracy_window) {
        ones=ones+ 1
      }
      
    }
    cat("\n Completed (no. trials): ",N)
    p <- ones / averagingN
    my_ps<-c(my_ps,p)
  }
  return(my_ps)
}
