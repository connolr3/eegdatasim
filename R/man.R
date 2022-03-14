
#See: https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator

library(R.matlab)
library(RandomFields)
library(optimr)


meanpower <- readMat("meanpower.mat")$meanpower

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
noise <- function(frames, epochs, srate) {
  sumsig = 50#number of sinusoids from which each simulated signal is composed of
  #signal <- matrix(0,1,epochs*frames)
  signal <- replicate(epochs * frames, 0)
  for (trial in 1:epochs) {
    freq = 0
    signal_range <- c(((trial - 1) * frames + 1):(trial * frames))
    for (i in 1:sumsig) {
      freq = freq + (4 * runif(1, 0, 1))
      freqamp = meanpower[min (ceiling(freq), 125)] / meanpower[1]#meanpower created up above from read-in file
      phase <- runif(1, 0, 1) * 2 * pi
      signal[signal_range] <-
        signal[signal_range] + sin(c(1:frames) / srate * 2 * pi * freq + phase) * freqamp

    }
  }
  signal
  return(signal)
}




#' @export
peak <-
  function(frames,
           epochs,
           srate,
           peakfr,
           position = frames / 2,
           tjitter = 0,
           wave = NULL) {
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





#reshape approach

#' @export
signal_averaging1<-function(data,frames,epochs){
  m1 <- matrix(data, ncol=frames, byrow=TRUE)
  d1 <- as.data.frame(m1, stringsAsFactors=FALSE)
  summations<-colSums(d1)
  average_signal<-summations/epochs
}

#approach by summation

#' @export
signal_averaging2<-function(data,frames,epochs){
  average_signal<-rep(0,frames)
  range<- c(1:frames)
  for (i in c(1:epochs)) {
    average_signal<-average_signal+data[range]
    range<-range+frames
  }
  return(average_signal/epochs)
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
  sig_hat<-sqrt(var(reg_data))
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
  buffer_range<-lo:hi
  reg_data<- data[-buffer_range]
  sig_hat<-sqrt(var(reg_data))
  normhat<-mean(reg_data)
  return(c(sig_hat,normhat))
}

#' @export
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



#plots signal with erp highlighted in red
plot_erp<-function(signal,erp_range){
  x_values<-c(1:length(signal))
  point_colour<-replicate(length(signal),"black")
  point_colour[erp_range]<-"red"
  mydf<-data.frame(x_values,signal,point_colour)
  plot(x_values,signal,col=point_colour,main="Signal with ERP identified")
}



min_SSE<-function(data,par,srate,peakcenter){
  #  pred_values<-(1+cos(2*pi*par[1]*(data["x"]-par[2])))/2
  pred_values<-cos(((data["x"]-peakcenter)*2*pi*par[1])/srate)
  errs<-data["y"]-par[2]*pred_values
  sum((errs)^2)
}

#' @export
optimise_ERP<-function(dat,startingfrequency=0,starting_amp=1,mysr,pkcntr){
  result <- optim(par = c(startingfrequency, starting_amp), fn = min_SSE, data = dat,srate=mysr,peakcenter=pkcntr)
  return(result)
}






