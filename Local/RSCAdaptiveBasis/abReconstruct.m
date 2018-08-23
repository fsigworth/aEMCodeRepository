function vol=abReconstruct(projs,ctfs,angles)
% function vol=abReconstruct(projs,ctfs,angles)
% projs: a stack of class-mean images n x n x nim
% ctfs: n x n x nim
% angles: nim x 3 array of rsc angles in degrees

data_axes = fredAnglesToAB(angles);

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
nest_plot_flag    = -1;          % Plot intermediate slices of iterations of Nesterov loop
%   -1 means suppress all output.
basis             = 's';        % 's' for swt16 frame, 'w' for regular wavelet basis
thr_type          = 0;          % 0 Soft Threshold, 1 Hard Threshold
W                 = 'coif3';    % Wavelet family
W_lev             = 2;          % Wavelet tree depth. Note: if basis = 's', then W_lev must = 2.

% Nesterov Algorithm Limits (can play with these if things don't converge)
nesterov_iter_lim = 15;
nesterov_stop_lim = 0.08;

% STEP SIZE for Nesterov 
mu                = 1e-4;       % Note: this will depend on your dataset

% Calculate some numbers
slice = ceil(size(projs,1)/2);
maskR = slice-1;



%% Run Adaptive Basis Reconstruction using Nesterov's Method
vol = reconstruct_by_nesterov_w_ctf(...
 		double(projs), data_axes, double(ctfs), maskR,...
 		basis, ...      
 		W, W_lev,thr_type,...
 		mu, nesterov_iter_lim, nesterov_stop_lim, fig_count, slice, nest_plot_flag);
