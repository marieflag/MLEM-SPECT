/* Description */
MLEM reconstruction for Dynamic Cardiac SPECT system developped in MGH&HMS. 
Bläckberg, L., Sajedi, S., Anderson, O.A., Feng, Y., El Fakhri, G., Furenlid, L. and Sabet, H., 2020, October. Dynamic Cardiac SPECT for diagnostic and theranostics applications: latest results. In 2020 IEEE nuclear science symposium and medical imaging conference (NSS/MIC) (pp. 1-3). IEEE.

DC-SPECT is a cardiac organ specific SPECT system that composed of multiple collimated heads, in this version the total number of heads is 80. This C++ code is parallelized with OpenMP for CPUs. The input data contains a MATLAB script describing the necessary parameters for setting up the reconstruction, a binary file containing projection on the 80 heads with a Derenzo source simulated with Gate V9.0, as well as a binary file containing simulated system matrix provided by collaborators in University of Arizona, using the GPU-accelrated Monte Carlo simulation software UAMC.
Kupinski, M.A., Ruiz-Gonzalez, M., Richards, R.G., May, M., Doty, K., King, M., Kuo, P. and Furenlid, L.R., 2022, November. Real-Time Modeling of the AdaptiSPECT-C Brain Imaging System for Hardware Evaluation, Acquisition Software Testing, and Adaptation-Rule Development. In 2022 IEEE Nuclear Science Symposium and Medical Imaging Conference (NSS/MIC) (pp. 1-2). IEEE.

Before execution, please download the simulated system matrix via the following link:
https://www.dropbox.com/scl/fi/qt45whw1j0u3de2wrys46/psf100.binary?rlkey=mnfwxlg4xgdc5wnfiq2zvd1iw&e=1&st=ysk4p2bd&dl=0

/* Getting started */
//Dependencies cmake_minimum_required(VERSION 2.8) OpenMP
//Installing and excutating
cmake CMakeList.txt 
make 
./SPECT SPECT_head_PSF_derenzo.m 

By executing the code, users will obtain binary files containing raw images after each iteration of MLEM, as well as a text file containing events used in the reconstruction.
Details regarding to the parameters setting can be found in the matlab script SPECT_head_PSF_derenzo.m
