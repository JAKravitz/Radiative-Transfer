# Radiative-Transfer
A repository and working directory for processing code, data, and documents for radiative transfer modeling of synthetic spectral aquatic datasets for machine learning development.  

## MODTRAN 
* RunBatch_v#.m is main processing script which calls and begins batch modeling of hydrolight output to TOA radiances
* Before RunBatch, FindHyCases.m be pointing to directory of hydrolight output and run separately to create HyCaseSetup script which compiles all batch hydrolight runs for Modtran
* ModCaseSetup.m should be previously setup for required atmospheric parameters
* ExecuteCase.m should also be pointing to hydrolight output directory, and performs the modtran heavy lifting. MATLAB wrapper for modtran requried.
* Other scripts used for resolving spectra to sensor SRFs and including sensitivity information  
