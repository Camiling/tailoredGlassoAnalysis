# tailoredGlassoAnalysis
Analysis of simulated and real data using the tailored graphical lasso. 

Code for performing the simulations is given in R/simulation_study.R, using the functions given in R/simulation_functions.R. Code for performing the extended simulation study is given in R/simulation_study_extended.R, using the functions given in R/simulation_functions.R. Note that for the extended simulation study, the package located at R/espace must be installed. 

Code for downloading and analysing the Oslo 2 data is given in R/analyse_Oslo2.R, and code for analysing the TCGA BRCA data is given in R/analyse_TCGA.R. Annotation data used in the analyses is given in annotation_data/. For each data set, a list of all edges in the tailored graphical lasso graph, as well as a list of the edges present in the tailored graphical lasso graph but not the weighted graphical lasso graph, is given in Edge_lists/. 

