#See: https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
usethis::use_package("R.matlab")

#' @export
plot_signal <- function(data, mark = NULL, title = "EEG Signal") {
  plot(data, main = title, type = "l")
  # this if statement and corresponding code is applicable for phase resetting theory only
  if (is.null(mark) == FALSE) {
    abline(v = mark, col = "red")
    text(mark + 5, 0.25, paste(mark), srt = 90, col = "red")
  }
}


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

#' @export
noise <- function(frames, epochs, srate, meanpower = NULL) {
  if (frames < 0) stop("frames cannot be less than 0")
  if (srate < 0) stop("srate cannot be less than 0")
  if (epochs < 0) stop("epochs cannot be less than 0")
  if( is.null(meanpower)){
  #  path<-file.path("R","meanpower.mat")
  #  meanpower <- as.vector(R.matlab::readMat(path)$meanpower)
    meanpower<-c(0.001512436 ,0.0008492235 ,0.0006087524 ,0.0004898673 ,0.0004235233
                 ,0.0003916184 ,0.0003627388 ,0.0003405692 ,0.0003431604 ,0.0003972602
                 ,0.0004257465 ,0.0003469311 ,0.0002942894 ,0.0002556115 ,0.0002355836
                 ,0.0002227425 ,0.0002139522 ,0.0002077288 ,0.0002048202 ,0.0002017179
                 ,0.0002007883 ,0.0001979438 ,0.0001955506 ,0.0001913121 ,0.0001866472
                 ,0.0001830944 ,0.0001782477 ,0.0001759883 ,0.0001716569 ,0.0001692953
                 ,0.0001641865 ,0.0001611862 ,0.0001595688 ,0.0001574246 ,0.0001553133
                 ,0.0001546942 ,0.0001521656 ,0.0001495667 ,0.0001471883 ,0.0001465332
                 ,0.0001451215 ,0.0001449061 ,0.0001434163 ,0.0001444918 ,0.0001419088
                 ,0.0001434826 ,0.0001430495 ,0.0001411656 ,0.0001425278 ,0.0001405317
                 ,0.0001406423 ,0.0001381596 ,0.0001383266 ,0.00013682 ,0.0001370514
                 ,0.0001373267 ,0.0001371059 ,0.000135919 ,0.0001367891 ,0.0001448403
                 ,0.0001355708 ,0.0001367028 ,0.0001360831 ,0.0001333923 ,0.0001329323
                 ,0.0001326268 ,0.0001309518 ,0.0001319445 ,0.0001302691 ,0.0001297129
                 ,0.0001288177 ,0.0001280689 ,0.000129027 ,0.0001311451 ,0.0001288131
                 ,0.0001293273 ,0.0001274819 ,0.000126366 ,0.0001278305 ,0.0001259032
                 ,0.0001262501 ,0.0001252715 ,0.0001256308 ,0.0001269599 ,0.0001253855
                 ,0.0001258603 ,0.0001232594 ,0.0001247888 ,0.0001235551 ,0.0001233012
                 ,0.0001229693 ,0.0001232385 ,0.0001215312 ,0.0001221649 ,0.0001226418
                 ,0.000121279 ,0.0001224354 ,0.0001222475 ,0.0001204993 ,0.0001210175
                 ,0.0001219137 ,0.0001208042 ,0.0001191548 ,0.0001192165 ,0.0001196114
                 ,0.0001190726 ,0.0001185056 ,0.0001184628 ,0.0001203116 ,0.0001189227
                 ,0.0001179185 ,0.0001202338 ,0.0001197265 ,0.0001200837 ,0.0001174336
                 ,0.0001177806 ,0.0001195453 ,0.0001210516 ,0.0001182723 ,0.0001175911
                 ,0.0001172347 ,0.0001172174 ,0.0001188854 ,0.0001175014 ,0.0001183316
    )
  } 
  signals<-matrix(0,frames,epochs)
  signals<-apply(signals,2,fill_signals,meanpower=meanpower,frames=frames,epochs=epochs,srate=srate)
  return(as.vector(signals))
}

#' @export
peak <- function(frames, epochs,srate,peakfr,position = frames / 2,tjitter = 0,wave = NULL) {
  if (frames < 0) stop("frames cannot be less than 0")
  if (srate < 0) stop("srate cannot be less than 0")
  if (epochs < 0) stop("epochs cannot be less than 0")
  mypeak <- c(1, 5, 4, 9, 0)
  signal <- replicate(epochs * frames, 0)
  for (trial in 1:epochs) {
    pos = position + round(runif(1, 0, 1) * tjitter)
      for (i in 1:frames) {
        phase = (i - pos) / srate * 2 * pi * peakfr
        if ((is.null(wave) == FALSE) ||(phase < pi / 2 && phase > -pi / 2)) {
          #if wave | (phase < pi/2 & phase > -pi/2)
          signal[(trial - 1) * frames + i] = cos(phase)
        }
      }
    }
    return(signal)
  }


#' @export
phasereset <- function (frames,epochs,srate,minfr,maxfr,position = frames / 2,tjitter = 0) 
  {
    signal <- replicate(epochs * frames, 0)#generating empty wave
    for (trial in 1:epochs) {
      wavefr = runif(1, 0, 1) * (maxfr - minfr) + minfr
      initphase = runif(1, 0, 1) * 2 * pi
      pos = position + round(runif(1, 0, 1) * tjitter)
      for (i in 1:frames) {
        if (i < pos) phase = i / srate * 2 * pi * wavefr + initphase
        else phase = (i - pos) / srate * 2 * pi * wavefr
        signal[(trial - 1) * frames + i] = sin(phase)
      }
    }
    return(signal)
  }


#' @export
Makinen <- function(frames, epochs, srate, position = frames / 2) {
  signal <- replicate(epochs * frames, 0)#generating empty wave
  for (i in 1:4) {#repeat for 4 sinusoids
    new_signal <- phasereset(frames, epochs, srate, 4, 16, position)
    signal = signal + new_signal#add new sin to summation
  }
  return(signal)
}


#' @export
Makinen1a <- function() {
  trials = 30
  mysignal = Makinen (400, trials, 1000, 175)#30 trials
  my_new_signal <- matrix(mysignal, nrow = 400, ncol = 30)#reshape into matrix
  my_df <-as.data.frame(t(my_new_signal))#convert -> df to use a matlab plot

  par(mfrow = c(3, 1))
  matplot(t(my_df), type = "l", ylab = "EEG")
  signalmean <- lapply(my_df[1:400], FUN = mean)
  plot( c(1:400),signalmean,type = "l",col = "blue",xlab = "",ylab = "ERP")
  variance <- lapply(my_df[1:400], FUN = var)
  plot(c(1:400),variance,type = "l",col = "blue",xlab = "")
  par(mfrow = c(1, 1))#reset par to normal
}

#' @export
signal_averaging<-function(data,frames,epochs){
  a <- matrix( data, nrow=frames, ncol=epochs )
  return(rowMeans(a))
}

#' @export
estimate_amplitude<-function(averaged_signal){
  return(max(averaged_signal))
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
#' @export
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
  result <- optim(par = 1, fn = min_SSE, gr=gr_min_SSE, x=x, y=y,srate=mysr,peakcenter=pkcntr, method="BFGS")
  z <- cos(((x-pkcntr)*2*pi*result$par[1])/mysr )
  alphaest <- sum(z * y) / sum(z*z)
  result$par <- c((result$par),alphaest)
  return(result)
}

#' @export
power_determination <-function(accuracy_window,freq,amp,frames,srate,maxtrial = 40,averagingN=100)
{
  set.seed("123")
  my_frequency_ps <- numeric(length=maxtrial)
  my_amplitude_ps <- numeric(length=maxtrial)
  meanpower <- c(0.001512436 ,0.0008492235 ,0.0006087524 ,0.0004898673 ,0.0004235233
                 ,0.0003916184 ,0.0003627388 ,0.0003405692 ,0.0003431604 ,0.0003972602
                 ,0.0004257465 ,0.0003469311 ,0.0002942894 ,0.0002556115 ,0.0002355836
                 ,0.0002227425 ,0.0002139522 ,0.0002077288 ,0.0002048202 ,0.0002017179
                 ,0.0002007883 ,0.0001979438 ,0.0001955506 ,0.0001913121 ,0.0001866472
                 ,0.0001830944 ,0.0001782477 ,0.0001759883 ,0.0001716569 ,0.0001692953
                 ,0.0001641865 ,0.0001611862 ,0.0001595688 ,0.0001574246 ,0.0001553133
                 ,0.0001546942 ,0.0001521656 ,0.0001495667 ,0.0001471883 ,0.0001465332
                 ,0.0001451215 ,0.0001449061 ,0.0001434163 ,0.0001444918 ,0.0001419088
                 ,0.0001434826 ,0.0001430495 ,0.0001411656 ,0.0001425278 ,0.0001405317
                 ,0.0001406423 ,0.0001381596 ,0.0001383266 ,0.00013682 ,0.0001370514
                 ,0.0001373267 ,0.0001371059 ,0.000135919 ,0.0001367891 ,0.0001448403
                 ,0.0001355708 ,0.0001367028 ,0.0001360831 ,0.0001333923 ,0.0001329323
                 ,0.0001326268 ,0.0001309518 ,0.0001319445 ,0.0001302691 ,0.0001297129
                 ,0.0001288177 ,0.0001280689 ,0.000129027 ,0.0001311451 ,0.0001288131
                 ,0.0001293273 ,0.0001274819 ,0.000126366 ,0.0001278305 ,0.0001259032
                 ,0.0001262501 ,0.0001252715 ,0.0001256308 ,0.0001269599 ,0.0001253855
                 ,0.0001258603 ,0.0001232594 ,0.0001247888 ,0.0001235551 ,0.0001233012
                 ,0.0001229693 ,0.0001232385 ,0.0001215312 ,0.0001221649 ,0.0001226418
                 ,0.000121279 ,0.0001224354 ,0.0001222475 ,0.0001204993 ,0.0001210175
                 ,0.0001219137 ,0.0001208042 ,0.0001191548 ,0.0001192165 ,0.0001196114
                 ,0.0001190726 ,0.0001185056 ,0.0001184628 ,0.0001203116 ,0.0001189227
                 ,0.0001179185 ,0.0001202338 ,0.0001197265 ,0.0001200837 ,0.0001174336
                 ,0.0001177806 ,0.0001195453 ,0.0001210516 ,0.0001182723 ,0.0001175911
                 ,0.0001172347 ,0.0001172174 ,0.0001188854 ,0.0001175014 ,0.0001183316
  )
  for (N in 1:maxtrial)
  {
    freq_ones <- 0
    amp_ones <- 0
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
      if (abs(freq - pars$par[1]) <= freq * accuracy_window) 
        freq_ones=freq_ones+ 1
      if (abs(amp - pars$par[2]) <= amp * accuracy_window) 
        amp_ones=amp_ones+ 1
    }
    
    cat("\n Completed (no. trials): ",N)
    freq_p <- freq_ones / averagingN
    amp_p <- amp_ones / averagingN
    my_frequency_ps[N]<-freq_p
    my_amplitude_ps[N]<-amp_p
  }
  results<-data.frame(my_frequency_ps,my_amplitude_ps)
  return((results))
}
