#run devtools::test()

testthat::test_that("errors", {
  mynoise<-noise(200,5,250)
  mypeak<-peak(200,5,250,7)
  mysignal<-mynoise+mypeak
  my_a_signal<-signal.averaging(mysignal,200,5)
  my_p<-power.determination(0.1,7,5,200,250,5,5)
  my_hats<-est.sig.hat(my_a_signal,100,0.2)
  my_stan_data<-(my_a_signal-my_hats[1])/my_hats[2]
  my_erp_range<-find.ERP.range(my_stan_data,2)
  my_estimations<-optimise.ERP(my_erp_range,my_a_signal[my_erp_range],250,100)
  
  #test error messages
  testthat::expect_error(
    noise(-200,5,250),
    "frames cannot be less than 0"
  )
  testthat::expect_error(
    noise(200,-5,250),
    "epochs cannot be less than 0"
  )
  testthat::expect_error(
    noise(200,5,-250),
    "srate cannot be less than 0"
  )
  testthat::expect_error(
    peak(-200,5,250),
    "frames cannot be less than 0"
  )
  testthat::expect_error(
    peak(200,-5,250),
    "epochs cannot be less than 0"
  )
  testthat::expect_error(
    peak(200,5,-250),
    "srate cannot be less than 0"
  )
  testthat::expect_error(
    peak(200,5,250,7,-4),
    "position must be in range 0:frames"
  )
  testthat::expect_error(
    estimate.amplitude("wrong type"),
    "averaged_signal must be in vector form"
  )
  testthat::expect_error(
    est.sig.hat(my_a_signal,100,0.55),
    "buffer_pc is too large. Choose a value <0.45"
  )
  testthat::expect_error(
    optimise.ERP(my_erp_range,c(my_erp_range,22),250,100),
    "lengths of x and y must equal"
  )
  false_x <- my_erp_range[!my_erp_range %in% 100]
  testthat::expect_error(
    optimise.ERP(c(false_x,22),my_stan_data[my_erp_range],250,100),
    "x must contain peak centre"
  )
  testthat::expect_error(
   power.determination(0.1,-5,5,200,250,5,1),
    "freq cannot be less than 0"
  )
  
  
  #test function return types and lengths
  testthat::expect_vector(mynoise,5*200)#noise return type +length
  testthat::expect_true(is.numeric(mynoise))
  
  testthat::expect_vector(mypeak,5*200)#peak return type +length
  testthat::expect_vector(mysignal,5*200)
  testthat::expect_true(is.numeric(mypeak))
  
  testthat::expect_vector(my_hats,2)#test est_sig_hat return type +length
  testthat::expect_vector(my_a_signal,200)#signal_averagingreturn type +length
  testthat::expect_true(is.numeric(my_a_signal))
  
  testthat::expect_length(my_estimations,5)
  testthat::expect_true(is.list(my_p))
  testthat::expect_length(my_p,2)
  testthat::expect_true(is.data.frame(my_p))
  
  #test screen output of print functions
  testthat::expect_output(power.determination(0.1,7,5,200,250,5,5))
})
