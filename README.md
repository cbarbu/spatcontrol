spatcontrol
===========

This package aims at exploring the determinants of spatial data. It focuses on complex spatial patterns: influence of known barriers on the dispersion/spatial autocorrelation, multiple spatial scale, etc... 

Functions are in extrapol\_field.R. Example of the use of the main functionalities are given the in the example\_\* files:
- generate are spatially autocorrelated data (see gen.map, zgen)
- generate structured autocorrelograms to examine the impact of known barriers
- fit a gaussian field accounting for known barriers, cofactors and observers quality in the fit.

Other utilities can be found in extrapol\_field.R, organized into chapters:
- General purpose functions
- Functions specific to the autocorrelation analysis
- Map generation
- Functions specific to the GMRF

The spatcontrol package is under development at the github repository: 
https://github.com/cbarbu/spatcontrol 

Participation are welcome through forking and pull request in GitHub. The aim is to contribute this work and others at R CRAN package. 
