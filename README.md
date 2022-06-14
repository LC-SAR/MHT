# MHT
MHT (Multiple Hypothesis Testing) is a means to identify the most probable mathematical model, and designed for InSAR (deformation) time series modeling. For each InSAR measurement point, the null hypothesis of ‘steady-state’ motion is considered as default, which is tested against a multitude of potential temporal models, built based on a library of canonical functions. If the null hypothesis is sustained, there is no (significant) anomaly in the data. If the null hypothesis is rejected, we test the entire library of potential alternative models with different physically realistic parameters against the null hypothesis using the B-method of testing. Finally, using test-ratios, we select the most likely model for each InSAR measurement point, update the quality description of the estimates, while avoiding overfitting. 

The MHT scripts are developed upon Matlab, and can be easily converted to Python format.  Note that this is the research code provided to you 'as it' with no warranties of correctness. Use at your own risk.


Citation:
Chang, L., Hanssen, R.F., 2015. A probabilistic approach for InSAR time-series postprocessing. IEEE Transactions on Geoscience and Remote Sensing. 54 (1), 421–430.
