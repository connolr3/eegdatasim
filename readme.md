EEGDataSim 
R version of https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator 
with additional functionality
__________________
DEPENDENCIES

The following packages are required for eegdatasim;
+ install.packages("R.matlab")
+ install.packages("devtools")   - required for testing only
__________________
INSTALLATION

To install this package run:
+ devtools::install_github("connolr3/eegdatasim")
__________________
SHINY APP

There is also a shiny app that acts tutorial to this packages. Please see it hosted on: https://connolr3.shinyapps.io/EEGDataSim/
The repo for this Shiny App is: https://github.com/connolr3/eegdatasim_shiny
__________________
TESTING

Test file (testthat.R) is included in tests folder. Navigate to test file in RStudio, ensure working directory is set to tests folder and run;
+ devtools::test()
