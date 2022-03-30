
testthat::test_that("errors", {
  #testnoise
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
 # testthat::expect_length(noise(200,5,250),5*200)
  #testpeak
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
  testthat::expect_length(peak(200,5,250,7),5*200)
  
  
  
})