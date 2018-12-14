function [spec, shot, gaussians]=NoiseModel3(f,p)
% Sum of Gaussians with p=amplitudes.
% A matrix nf x np of individual Gaussians is returned if desired. 
sz=size(f);
f=f(:);
p=p(:);
f0=.01; % lowest frequency that concerns us
np=numel(p)-2;
k=-0.5./(f0*2.^(0:np-1)).^2;
gaussians=exp(f*k); % nf x np matrix
spec=reshape(gaussians*p,sz);
shot=p(np)*(1-p(np+1)*;
