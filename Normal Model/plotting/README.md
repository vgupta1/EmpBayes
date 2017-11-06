## Overview

This folder contains various R-scripts mostly used to create pretty plots of simulation results.  
For the most part, actual simulation and competition use the *.jl files in the parent directory.  
Hence, these files are not strictly necessary to reproduce experimental outcomes.  

The most critical files in this folder are:
 - CLTPlots.R  - Used to generate the various "Robustness to Nonnormality" plots in the paper.
 - LOO_ThmPlot.R - Used to generate the plot for Theorem 4.1.  
 - tauDependence.R - Used to generate Figure 2a "Benefits over Estimate then Optimize"
 - plot3Part.R - Used to generate Figure 2b and 3: "Effects of Shrinkage on Solution to Example 3.1"
 - portfolioPlots.R - Used to generate all plots pertaining to the POAP case study. 
 - simSummaries.R - Used to generate the summary histograms for the simulation method of the POAP case study.  
 
 Other files pertain to plots that were not ultimately included in the final manuscript.  
