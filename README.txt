
To run the reproduction scripts first install the R package hyperlatent either directly from the accompanying tarball 

hyperlatent_1.0.tar.gz

or directly from github using devtools 

R> library(devtools)
R> install_github("grosed/hyperlatent")


The reproduction scripts and associated data are contained in the accompanying zip file

 reproduction.zip
 
 Unzip this file to run the replication scripts as outlined in the following :
 

1) 'reproduction/assessing_model_fit' contains code to implement the study in Section 6.1. The study is run by first executing 'modelfit.R' and plots are created by executing 'plot_modelfit.R'.

2) 'reproduction/misspecification' contains code to implement the study in Supplement K. The study is run by first executing 'misspecsims.R' using the bash script 'run_misspec.sh', and plots are created by executing 'plot_misspecout.R'.

3) 'reprodution/model_depth_comparisons' contains code to implement the study in Supplement H. The study is run by first executing 'run_modelcomp.R' and plots are created by executing 'plot_modcomp.R'

4) 'reproduction/scalability' contains code to implement the study in Section 6.2. The study is run by executing 'scalability_study.R' and plots are created by executing 'plot_scalability.R'

5) 'reproduction/real_data_examples' contains code to implement the data examples in Sections 7.1 and 7.2. A description of workflow for each example is given below. 

a) For the grocery example in Section 7.1, the data are loaded from the arules R package. The script 'reproduction/real-data-examples/data/process_grocery_data.R' is first run to load in the data and create a subsample as described in the main article. The subsampled data are saved in the folder 'reproduction/real-data-examples/data' and, once created, model fitting is then run using the 'fit_to_grocery.R' script. 

b) For the DBLP data example in Section 7.2, the script '.reproduction/real-data-examples/data/process_dblp.R' has been run, the subsampled data will be saved in the folder 
'reproduction/real-data-examples/data'. Then, model fitting is implemented using the 'fit_to_dblp.R' script.

Finally, plots for both data examples are created by executing 'plot_data_predictives.R'.
