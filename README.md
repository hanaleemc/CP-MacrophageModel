# CP-MacrophageModel
Included is my code for generating a Macrophage Model from RECON 3D and the Sanity Tests to assess basic model performance. 

Included in this reposity are the following files: 
Msample_analysis_full.R
macrophage_seedmodel
sanity_tests.m

Msample_analysis_full.R contains all the data cleaning/normalization code used to extract gene expression data through Agilent or Affy analysis. This code transforms the expression data into a series of 0s(inactive) and 1s(active) this information forms the basis of my construction of a macrophage model. 

Macrophage_seedmodel is a MATLAB file which generates a macrophage model from RECON3D, our most updated model of human metabolism. This script uses the expression data obtained in Msample_analysis_full.R to select only active genes from RECON3D as the basis of the context-specific model of macrophage cells.

Sanity_tests.m is a MATLAB file which goes through the Heirendt et al. (2019) COBRA toolbox ‘Testing Basic properties of a metabolic model” in order to ensure that the model is meeting basic quality checks. This is a series of 12 MATLAB tests, which assesses different functions of the model, to ensure it is meeting basic standards.
