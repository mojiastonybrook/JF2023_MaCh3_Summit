# JF2023_MaCh3_Summit
This documetation describes the general work flow of the joint analysis of T2K and SK neutrino experiments in 2023 with a Bayesian fitter, known as MaCh3, on the compute cluster Summit.
## Install MaCh3 on Summit
MaCh3 repository locates at:`https://github.com/t2k-software/MaCh3`. You need to be added to the member list to access the page. 
The branch used in joint fit 2023 is: `DBarrow_JointFit_CorrelatedFDDetSysts_MinorUpdates`

A reference guide on how to install MaCh3 on the Summit cluster: `https://github.com/weishi10141993/MaCh3/blob/DBarrow_JointFit/MaCh3_Summit.md#mach3-installation-on-summit`
## Massive production of MCMC on Summit
Based on the compute resource architecture and job management system adopted by Summit, a supportive package of scripts is used to prepare the necessary files for generating MCMC samples in parallel: `https://github.com/mojiastonybrook/summit_mach3_jobscripts`.

In the template `https://github.com/mojiastonybrook/summit_mach3_jobscripts/blob/master/run_mach3_TEMPLATE.sh#L7C9-L7C9` specify the directory to the one where MaChs is installed.

The configuration files in the directory of `SampleConfigs` need to be checked before using to make sure they are consistent with the tuning setup of interest for MaCh3, especially the atmospheric configuration card `https://github.com/mojiastonybrook/summit_mach3_jobscripts/blob/master/SampleConfigs/AtmConfig.cfg`.

The executable in MaCh3 to be run is specifiied by the command:`https://github.com/mojiastonybrook/summit_mach3_jobscripts/blob/master/prepare_job_scripts.py#L155`. 
The Asimov fit is done by `JointAtmFit`; the data fit is done by `JointAtmFit_Data`. 

After running the python script, the configs and directories for the chain production would be generated, along with the batch and run scripts to excuate the iterative prodcutions. 
### Instant rough check on the status of chain production
`https://github.com/mojiastonybrook/JF2023_MaCh3_Summit/blob/master/production_check_mod/check_JointAtmFit_outputs.py`

This script will examine the existence of the output root files from each chain. 
## Reduction and merge of the chain

## Reweighting with reactor constraints

## Making contour plots with the reduced chain




