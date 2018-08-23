function w = get_wavelet_grad_with_ctf(s_wdec,maskR,proj,coord_axes,ctfs)
n_proj = size(proj,3);

% Apply forward wavelet transform
s = waverec3(s_wdec);
g = zeros(size(s));

% Grab fields from wavelet structure
W_lev = s_wdec.level;
W     = s_wdec.filters;

% Check MATLABPOOL state and initiate 4-core processing
if(matlabpool('size')==0)
  matlabpool local
end

parfor i=1:n_proj,
  g = g + ...
    get_grad_component(s, maskR, ...
    squeeze(proj(:,:,i)),...
    coord_axes(:,i),...
    squeeze(ctfs(:,:,i)));
end

% Apply backward wavelet transform
g_wavelet_comp = wavedec3(g,W_lev,W,'mode','per');
w              = g_wavelet_comp.dec;

function g_comp = get_grad_component(s,maskR,class_mean,cur_axes,ctf)
% Forward project and apply CTF
img = mex_forward_project(s,cur_axes,maskR);
img = fftshift(fft2(img));
img = img.*ctf;
    
% Calculate difference
img = fftshift(fft2(class_mean)) - img;

% Apply CTF and backproject
img    = img.*ctf;
img    = real(ifft2(ifftshift(img)));
g_comp = mex_back_project(img,cur_axes,maskR);


