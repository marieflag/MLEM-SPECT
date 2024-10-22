% This script contains all the parameters needed for compiling a customized MLEM C++ code using UAMC simulated
% system matrix. Users can change parameters in: Input, Output, Recon volume, and Iterations.
% The RL-deconvolution and Total Variation regularization are turned off.
% Analytical calculation of the system matrix is not included in this version of MLEM.
%% Input

DATA_FILE = '/home/rpiladmin/CodeC++_multihead_Matt/SPECT_head_PSF_derenzo4.sample0.Head80.bin'; % Projection on detectors, 1x50000 double (25x25x80)
SIMU_FILE = '/home/rpiladmin/CodeC++_multihead_Matt/37tripsf150.binary'; % Simulated system matrix

%% Output

RESULTS_DIR = '/home/rpiladmin/Desktop/'; % Path where the results are saved

%% Recon volume 
% Using the same size as the simulation for the system matrix

VOLUME_DIMENSIONS = [15, 15, 15]; % in cm
VOXELS = [150, 150, 150]; 
VOLUME_CENTRE = [0, 0, 0];

%% Iterations

FIRST_ITERATION = 0; % If 0, starts from 0; if 10, starts from the 10th iteration
ITERATIONS = 400;
LAST_ITERATION = FIRST_ITERATION + ITERATIONS;

%% Total variation regularization

ALPHA_TV = 0; % TV regularization under construction, recommend setting to 0 to turn off the TV
TV_ITERATION = 20; % Iterations of TV

%% Camera properties 

NB_CAMERAS = 80;
DET_VOXELS = [25, 25, 1]; % 25x25 pixelated crystal per head
PINHOLE_SIZE = 0.22; % Unit cm, for analytic calculation of the system matrix, not used in the code
DET_SIZE = [5.6, 5.6, 2]; % Unit cm, not used in the code

%% The following four parameters should not be changed

DATA_TYPE = 'Gate';      
PSF = 'OFF'; % For analytic calculation of the system matrix, if 'ON' requires geometry description of the system
PSF_DECONV = 'OFF'; % RL-deconvolution, recommend setting to 0, if 'ON' filter files are required
VOXELS_FILTER = [11, 11, 11]; % Size of filter for RL-deconvolution, not used if PSF_DECONV is off
