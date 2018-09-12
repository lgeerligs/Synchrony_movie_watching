# Synchrony_movie_watching
This code can be used to reproduce the results that were reported in the Neurobiology of Aging paper “Age-related differences in information processing during movie watching”.  https://doi.org/10.1016/j.neurobiolaging.2018.07.025

This code reproduces the measures that were reported in the paper, performs the statistical analyses and produces most of the figures in the paper. It starts from the pre-processed and cleaned data, after all the steps described in the section ‘Functional data pre-processing’.


## Code
**compute_synchrony.m** computes the inter-subject synchrony for each subject. It also computes the correlation between synchrony and age, and performs the computation of Bayes factors. It also includes a number of control analyses, as described in the supplementary materials and the section “Sensory loss”. This script generates a number of files that are used as input files for the other scripts below. 

**relate_synchrony_to_cognition.m** relates synchrony to cognitive performance scores. 

**WSBM_community_detection.m** find communities of participants based on the timecourses within clusters.

**synchrony_over_time.m** computes time-varying synchrony using sliding window analyses. 

**relate_synchrony_to_connectivity.m** relates synchrony to functional connectivity. It also computes the functional connectivity summary scores and relates those to the white matter (MK) values.  

Helper scripts that are required in parts of the code are contained in a separate folder. The folder ‘external’ contains helper scripts that were written by others and that I either downloaded via MatlabCentral or found online somewhere else. The Bayesfactor toolbox included here is the Matlab version of the R-toolbox https://github.com/MicheleNuijten/BayesMed. I found this Matlab version online but it is no longer available, which is why I included it here for the sake of reproducibility. 
Other toolboxes used in these scripts are:
SPM12 - https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
Mediation Toolbox - https://github.com/canlab/MediationToolbox 
Canlab core toolbox - https://github.com/canlab/CanlabCore 


## Data 
The data files descibed below are available from 11633/di.dcc.DSC_2018.00013_064 

**residuals_subject_selection.mat** contains the a number of attributes of the participants, including their age, gender, head motion (max - rms_max and total - rms_tot) and subject IDs (CCID). It also contains the variable residuals, which contains the timecourses of each of 748 ROIs and each participant (in Craddock_ROIs_included_ordered.nii) 

**Modules.mat** contains the variable modcor, which gives the module numbers for each of the ROIs in the data. The network labels correspond to each of the module numbers of in the variable network_name. ROIs numbered with a 0 were not part of any of the modules. 

**Behavioral_data.mat** contains the cognitive performance scores for all participants

**WMdata.mat**, the variable WMdata contains the mean kurtosis values for a set of WM tracts from the JHU atlas (JHU-ICBM-tracts-maxprob-thr25-2mm.nii), it also contains the intracrantial volume (TIV) and an index of the stripes in the diffusion data (stripeind). Participants with values larger than 0.1 were excluded from the analysis. 

