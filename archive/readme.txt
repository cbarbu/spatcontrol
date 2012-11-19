July 26, 2012 - reverted to the simple log sampling as the "centered" log sampling was doing weird things and theoretically the simple log sampling should return correct values (only think needed on the sampling function is to set well the hasting term to have reciprocal values for theta and theta^star. Not symmetry of the sampling arround the starting value.

Aug 22, 2011 - Kc is reverted to a precision as expected through it's name, carefull 20110819-184901sub1PauTrue5Cof_NoNAInsp_eps_meanNorm_0.01_exp_Kc0.01 correspond to a very small variance.

# this folder is renamed as working, old version should be renamed, but working should keep its name
the main file for the fit is "full_sampler.r"
# can be run from R with "source()" and from bash by sec_launch.sh (in ~/bin/)

the file "pseudo_data_generation.r" can be used separatly for data generation or more generally for one time analysis of the dataset.

# the blockSize information is prepared in prep_impact_size_block.R
#

