%% Adaptive Basis Reconstruction %%

%% Initialize MATLAB and Set Parameters
close all
clear all
clc
warning('off','Images:imshow:magnificationMustBeFitForDockedFigure');
fig_count = 1;

% Check MATLABPOOL state and initiate multi-core processing
% Note: configuration 'local' should be properly configured for YOUR machine
% % if(matlabpool('size')==0)
% %      matlabpool local
% % end
% 
% Set path to utilities
% % path('utils',path);

% Nesterov Parameters (no reason to play with these)
nest_plot_flag    = 0;          % Plot intermediate slices of iterations of Nesterov loop
basis             = 's';        % 's' for swt16 frame, 'w' for regular wavelet basis
thr_type          = 0;          % 0 Soft Threshold, 1 Hard Threshold
W                 = 'coif3';    % Wavelet family
W_lev             = 2;          % Wavelet tree depth. Note: if basis = 's', then W_lev must = 2.

% Nesterov Algorithm Limits (can play with these if things don't converge)
nesterov_iter_lim = 15;
nesterov_stop_lim = 0.08;

% STEP SIZE for Nesterov 
mu                = 1e-4;       % Note: this will depend on your dataset

%% Load Example Inputs from Fred's Data Structure and Convert to Hemant's Notation
pa=ParsePath(which('adaptiveBasisReconstructionFRED'));
pa=[pa 'dataFromFred/'];
load([pa 'AlpsGoodies.mat']);

% Get the correct volume (for sanity checking)
correctVol = origVols;
clear origVols

% Get the images (or class means)
proj = double(imgsOrig);				% MUST BE DOUBLE, or HEMANT'S MEX FILES GO BONKERS
clear imgsOrig

% Use 'ones' image as CTFs for now
ctfs = ones(size(proj));

% Calculate some numbers
slice = ceil(size(proj,1)/2);
maskR = slice-1;

% Convert Fred's angles to Hemant's convention
data_axes = fredAnglesToAB(allSimAngles);
clear allSimAngles

%% Run Adaptive Basis Reconstruction using Nesterov's Method
tic
x_est = reconstruct_by_nesterov_w_ctf(...
 		proj, data_axes, ctfs, maskR,...
 		basis, ...      
 		W, W_lev,thr_type,...
 		mu, nesterov_iter_lim, nesterov_stop_lim, fig_count, slice, nest_plot_flag);
toc

%% Write output to a SPIDER volume
writeSPIDERfileALP('x_est.spi',x_est);
writeSPIDERfileALP('correctVol.spi',correctVol);
