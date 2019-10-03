# Convolutional neural network for seismic impedance inversion

Convolutional neural network for seismic impedance inversion

Paper reference:

Das, V., A. Pollack, U. Wollner, and T. Mukerji, 2019, Convolutional neural network for seismic impedance inversion: Geophysics, 84(6), 1-66.

(https://library.seg.org/doi/pdf/10.1190/geo2018-0838.1)

-- The folder ABC_uncertainty contains the work done for uncertainty quantification using Approximate Bayesian Computation. Jupyter notebook CNN_impedance_inversion-ABC.ipynb contains all the code related to computations. 

-- The folder Challenge_testing_Facies_variograms is one of the example folder containing all the input and output data for performing challenge testing using different values of the vertical range in the facies variograms. Jupyter notebook CNN_impedance_inversion.ipynb contains all the code related to computations. 

Generating the synthetic training dataset for the Challenge testing is done using the following MATLAB scripts and SGEMS (http://sgems.sourceforge.net/) software.

1. GetDifferentVariograms_Step0.m - Matlab script to get SGEMS script ready with different values of input variogram ranges 

2. reshapeImportedData_Step1.m - Matlab script to import results from SGEMS and reshape for further calculations 

3. RockPhysics_Step2.m - Matlab script to generate elastic properties from petrophysical properties

4. ForwardModel_Kennett_Step3.m - Matlab script to perform seismic forward modeling 

5. script_to_convert_data_Step4.m - Matlab script to convert all data and make it ready for use in python Jupyter notebook.

ForwardModel_Kennett_Step3.m uses a few scripts from SeisLab_10.0301 (https://www.mathworks.com/matlabcentral/fileexchange/15674-seislab-3-01)

-- The folder Volve_field_example contains code of using CNN for impedance inversion for Volve field data. Jupyter notebook Train_based_on_Volve_data.ipynb contains all the code related to computations. 



Note:

The dataset for Approximate Bayesian Computation is not included in the repository and can be obtained on request


