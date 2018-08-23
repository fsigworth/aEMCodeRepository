function [p T cc]=meMergeAligner2(m1,m2,cfilter,P,Psteps,niters,displayOn)
% function [p T]=meMergeAligner2(m1,m2,cfilter,P,Psteps,niters,displayOn)
% Optimize the alignment of two images using Simplex. The images are
% subjected to an affine transform (based on the parameters P
% through the function meMakeAffineMatrix) and then a translational
% cross-correlation is made via FFT, using the cfilter function (zero
% frequency is in the center) as a filter before inverse FFT.
% Given the transformation parameters P, vary them with initial steps
% given by Psteps (zero means a parameter not varied at all) and taking up to
% niters iterations.  If displayon=1, the cross-correlation is displayed as
% the fit is running.  The final parameter values and the final affine
% transformation matrix T are returned.  T(m2) matches m1.
% cfilter=filt;
% % cfilter=1;
% m1=SquareWindow(1024,64).*md1;
% subplot(2,2,2); imacs(m1);
% % r=radius(1024);
% 
% m2=md2;
% niters=10;
% displayOn=1;
% Psteps=Psteps1;


epsi=1e-5; % convergence criterion for Simplex
[ndx ndy]=size(m1);
nd=[ndx ndy];

nds=128;  % local align display size

fc1=fftn(m1).*cfilter;  % FT of reference
%%
done=0;
iters=0;
% P
p=Simplex('init',P,Psteps,Psteps>0);
T=meMakeAffineMatrix(p);
%%
% while iters<10 || ~(done || iters>niters)
while ~(done || iters>niters)
    m2s=AffineTransform(m2,T);
    ccd=fftshift(real(ifftn(fc1.*conj(fftn(m2s)))));
    [mx xs ys]=max2di(ccd);
    xsh=double(xs-(ndx/2+1))/ndx;
    ysh=double(ys-(ndy/2+1))/ndy;
    if displayOn && mod(iters,5)==0  % show the cc function
         ccalign=real(ifftn(fftn(ccd).*FourierShift(nd,-[xsh ysh].*nd)));
         cc=Crop(ccalign,nds);
        imacs(cc);  % show it with precise shifting
%          imacs(Crop(circshift(ccd,-round([xsh ysh].*nd)),nds));
        title(num2str([iters xsh*ndx ysh*ndy]));
        pstr=sprintf('%8.5f   ',p(1:4));
        xlabel(pstr);
        drawnow;
    end;
    %%
    p=Simplex(-mx);
    T=meMakeAffineMatrix(p);
    done=Simplex('converged',epsi*abs(mx));
    iters=iters+1;
end;
p=Simplex('centroid');
% update the translation parameters
p(5)=xsh;  % fractional shift of m2 to match m1
p(6)=ysh;
T=meMakeAffineMatrix(p);
