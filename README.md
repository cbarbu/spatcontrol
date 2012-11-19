spatcontrol
===========

This package aims at exploring the determinants of spatial data. It focuses on complex spatial patterns: influence of known barriers on the dispersion/spatial autocorrelation, multiple spatial scale, etc... 

Functions are in extrapol\_field.R. Example of the use of the main functionalities are given the in the example\_\* files:
- how to generate spatially autocorrelated data in example\_generation.R with gen.map
- how to compute structured autocorrelograms to examine the impact of known barriers on presence absence data in example\_structuredMI.R. 
- fit a gaussian field accounting for known barriers, cofactors and observers quality in the fit in exampled_fit_GMRF.R The parameters for the priors are configured in parameters\_extrapol.r

Other utility functions can be found in extrapol\_field.R, organized into chapters:
- General purpose functions (data management, sparse matrices handling, plotting)
- Functions specific to the structured autocorrelograms
- Map generation
- Functions specific to the GMRF

The spatcontrol package is under development at the github repository: 
https://github.com/cbarbu/spatcontrol 

Participation are welcome through forking and pull request in GitHub. The aim is to contribute this work and others at R CRAN package. 
