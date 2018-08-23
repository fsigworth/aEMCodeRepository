function s_est = reconstruct_by_nesterov_w_ctf(...
  proj, data_axes, ctfs, maskR,...
  basis, ...
  W, W_lev, thr_type,...
  mu, nesterov_iter_lim, nesterov_stop_lim, fig_count, slice, plot_flag)
%reconstruct_by_nesterov_w_ctf Adaptive Basis Reconstruction using Nesterov
%   s_est = reconstruct_by_nesterov_w_ctf(...
%       proj, data_axes, ctfs, maskR,...
%       basis, ...
%       W, W_lev, thr_type,...
%       mu, nesterov_iter_lim, nesterov_stop_lim, fig_count, slice, plot_flag)
%
%   proj: class mean images (stack)
%   data_axes: projection directions
%   ctfs: ctf images associated with each class mean image
%   maskR: radius to be used for spherical masking during reconstruction
%   basis: 's' for stationary frame, 'w' for wavelet basis
%   W: wavelet family name. We suggest 'coif3'
%   W_lev: wavelet level. We suggest 2
%   thr_type: threshold type. 0 for soft, 1 for hard thresholding
%   mu: step size for Nesterov
%   nesterov_iter_lim: limit for Nesterov iterations
%   nesterov_stop_lim: stopping threshold for Nesterov iterations
%   fig_count: number at which figures should begin
%   slice: slice through volume to be displayed
%   plot_flag: 0 for no images, 1 for slices at odd iterations.
%%%%%%%%%%     -1 for no printout on iterations

%   A.Kucukelbir 05-Jan-2011.
%   Last Revision: 15-Feb-2012.

%% MAD noise std and SURE threshold estimation

% Get the noise from the corners of the class means
n             = size(proj,1);
[dim1 dim2]   = ndgrid(-n/2:n/2-1); 
R             = sqrt(dim1.^2 + dim2.^2);
R             = (R>maskR);
R             = repmat(R,[1 1 size(proj,3)]);
noise         = proj(R);

% Calculate the MAD noist std and Lambda threshold
mad_std       = mad(noise,1)/0.6745;
lambda        = 3.0 * 2*sqrt(2) * mad_std * sqrt(2*log(n^3));
if plot_flag>=0
fprintf('reconstruct_by_nesterov_w_ctf: mad_std = %f\n',mad_std);
fprintf('reconstruct_by_nesterov_w_ctf: lambda  = %f\n',lambda);
fprintf('reconstruct_by_nesterov_w_ctf: lambda calculated as 3.0 * 2*sqrt(2) * mad_std * sqrt(2*log(n^3))\n');
end;

%% Initialize variables for reconstruction

% Create zero_struct (compatible with all bases)
tmp       = zeros(n,n,n);
zero_w    = wavedec3(tmp,W_lev,W,'mode','per');
zero_s    = swt16dec3(tmp,W_lev,W,'mode','per');

switch basis
  case {'w'}
    zero_struct            = zero_w;
  case {'s'}
    zero_struct            = zero_s;
    zero_struct.filterName = W;
end
clear tmp zero_w zero_s 

% Initialize parameters
xi      = zero_struct;
alpha   = zero_struct;
alpha_0 = alpha;
theta   = 0;
v       = zero_struct;
w       = zero_struct;
s_prev  = zeros(n,n,n);
eps_den = 1e-20;

%% Nesterov loop
for i = 1:nesterov_iter_lim
  
  % Calculate v
  v.dec = prox_regularizer( ...
    cellfun(@minus, alpha_0.dec, xi.dec,'Un',0),...
    theta*lambda, thr_type);

  % Calculate a
  a = (mu + sqrt(mu^2 + 4*mu*theta))/2;
  
  % Calculate w
  w.dec = cellfun( @(f) f./(theta+a),...
    cellfun(@plus,...
    cellfun(@(f) f.*theta, alpha.dec,'Un',0),...
    cellfun(@(f) f.*a,     v.dec,    'Un',0),...
    'Un',0),...
    'Un',0);
  
  % Calculate alpha
  switch basis
    case 'w'
      alpha.dec = prox_regularizer( ...
        cellfun(@plus, w.dec,...
        cellfun(@(f) f.*(0.5*mu), ...
        get_wavelet_grad_with_ctf(w,maskR,proj,data_axes,ctfs),...
        'Un',0),'Un',0),...
        0.5*mu*lambda, thr_type);
    case 's'
      alpha.dec = prox_regularizer( ...
        cellfun(@plus, w.dec,...
        cellfun(@(f) f.*(0.5*mu), ...
        get_stat_grad_with_ctf(w,maskR,proj,data_axes,ctfs),...
        'Un',0),'Un',0),...
        0.5*mu*lambda, thr_type);      
  end
  
  % Calculate xi
  switch basis
    case 'w'
      xi.dec = cellfun(@minus, xi.dec,...
        cellfun(@(f) f.*a, ...
        get_wavelet_grad_with_ctf(alpha,maskR,proj,data_axes,ctfs),...
        'Un',0),...
        'Un',0);
    case 's'
      xi.dec = cellfun(@minus, xi.dec,...
        cellfun(@(f) f.*a, ...
        get_stat_grad_with_ctf(alpha,maskR,proj,data_axes,ctfs),...
        'Un',0),...
        'Un',0);      
  end
  
  % Update theta
  theta = theta + a;
  
  % Calculate progress and breaking threshold
  switch basis
    case 'w'
      s_est = waverec3(alpha);
    case 's'
      s_est = swt16rec3(alpha);
  end
  
  % Display slices
  if plot_flag==1
    if mod(i,2)==1
      figure(fig_count)
      set(fig_count, 'WindowStyle', 'docked')
      imshow(s_est(:,:,slice),[])
      fig_count = fig_count + 1;
      drawnow
      pause(.5)
    end
  end
  
  % Calculate stopping threshold
  s_diff      = s_est - s_prev;
  s_diff_norm = norm(s_diff(:));
  s_norm      = norm(s_est(:));
  ds          = s_diff_norm/(s_norm + eps_den);
  
  if plot_flag>=0
      fprintf('reconstruct_by_nesterov_w_ctf: iter = %d ds=%f\n',i,ds);
  end;
  
  if(ds < nesterov_stop_lim)
    break
  end
  
  s_prev = s_est;
  
end


