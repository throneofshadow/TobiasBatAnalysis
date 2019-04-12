# TobiasBatAnalysis
Matlab code for bat call cutting and analysis

findcalls_v5_CoEd is the main pre-shift method for cutting calls.
findcalls_session_tobias is a function that also cuts preshift calls using Maimon's (maimonr) call cutting routine.
findcalls_tobias is the cutting function for post-shifted calls using Maimon's call cutting routine.
callData is the main analysis package. For more information on the function itself (non-modified function), check maimonr github page under callData repository. This version of callData has additional parameters added in from the generic version, including values for high/low pass shift values. This version is compatible with the workflow used in VocOperantTraining by Tobias Schmid. 

Jupyter notebooks are individual-level analysis done on selected calls from dataset. These include PCA/LDA analysis, as well as a random forest analysis for classification. Please ensure dependencies are in place before running. 
