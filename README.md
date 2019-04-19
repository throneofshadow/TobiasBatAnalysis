# TobiasBatAnalysis
Matlab code for bat call cutting and analysis

findcalls_v5_CoEd is the main pre-shift method for cutting calls. This function was written to be used in pre-sliced .mat files generated before the VocOperant System.
findcalls_session_tobias is a function (wrapper) that modifies VocOperant generated data into vocalizations using Maimon's (maimonr) call cutting routine. This is the function that you call to cut the calls. This function calls into findcalls_tobias, explained below.
findcalls_tobias is a unmodified version of Maimon's call cutting function. This is the function that actually cuts the calls.
callData is the main analysis package that generates a data object containing relevant parameters for each vocalization found using the workflow above. For more information on the function itself (non-modified function), check maimonr github page under callData repository. This version of callData has additional parameters added in from the generic version, including values for high/low pass shift values. This version is compatible with the workflow used in VocOperantTraining by Tobias Schmid. 

Jupyter notebooks are individual-level analysis done on selected calls from dataset. These include PCA/LDA analysis, as well as a random forest analysis for classification. Please ensure dependencies are in place before running. 
