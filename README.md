# Hybrid_Imputation_for_Missing_Values

## A repo for evaluating missing value imputation methods found in the R package imputeLCMD

*Miranda L. Gardner*    
*November 10, 2020*

## I. OVERVIEW
#### The purpose of this repository is to provide the database search results and R scripts utilized during the analysis and publication of MS data. This specific analysis was performed with two different types of DDA bottom-up mass spectrometry experiments, whole cell proteomics and immunoprecipitation. The public data was downloaded from PRIDE, searched with X!Tandem Fido on the OpenMS platform, analyzed with 6 types of missing value imputation methods and performance evaluated with statistical metrics following differential expression analysis.  
  
  
## II. INSTRUCTIONS FOR LAUNCHING BINDER AND RSTUDIO FOR DATA ANALYSIS

1. Click on the link below to launch the Binder Github Notebook. This environment is version-controlled (R 3.6.3) with the runtime.txt file found in this repository to ensure that all required packages used in this analysis are compatible.   
  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mlgred1205/Hybrid_Imputation_for_Missing_Values/HEAD)
  
2. Once Jupyter has opened, select a new RStudio session by clicking on the drop-down menu under the New tab at the right hand side of the screen.  
  
  
## III. DATA ANALYSIS  
1. Click on one of the three ".R" files to begin the data analysis and this will open the script in a window above the R console.  
  
2. To run the script, make sure the blinking cursor is in that window and hold control (or command for Mac) + A to select all. Once script is highlighted in blue, hit the run button (green arrow) above the script window.  
  
3. Note that this data analysis performs 25 consecutive imputations within a loop and there are quality control figures housed within this loop that are overwritten. In order to view these figures, you will need to select the lines of code within the loop and run those separately.  
  
  + MDA-MB-468 - select lines 335-536 and hit run  
  + EZH2toIgG - select lines 135-335 and hit run 
  + SUZ12toIgG - select lines 135-338 and hit run  
