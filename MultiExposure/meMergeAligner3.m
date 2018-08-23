function [p, T, cc]=meMergeAligner3(m1,m2,cfilter,P,Psteps,niters,displayOn)
% function [p T]=meMergeAligner3(m1,m2,cfilter,P,Psteps,niters,displayOn)
% Optimize the alignment of two images using Simplex. The images are
% subjected to an affine transform (based on the parameters P [1x6]
% through the function meMakeAffineMatrix) and then a translational
% cross-correlation is made via FFT, using the cfilter function (zero
% frequency is in the center) as a filter before inverse FFT.
% Given the transformation parameters P, vary them with initial steps
% given by Psteps (zero means a parameter not varied at all) and taking up to
% niters iterations.  %
%
% If P has two rows, then a brute-force search is made between P(1,:) and
% P(2,:) with steps Psteps.
%
% If displayon=1, the cross-correlation is displayed as
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
[ndx, ndy]=size(m1);
nd=[ndx ndy];

nds=round(nd/8);  % local align display size
mxv=-inf;

fc1=fftn(m1).*ifftshift(cfilter);  % FT of reference, times filter
%%
done=0;
iters=0;
%
if size(P,1)>1  % we are doing brute-force
    pv=P(1,:);  % initialization of optimized params
    activeDims=find(Psteps>0);
    nDims=numel(activeDims);
    if nDims<1
        error('No dimensions to search');
    end;
    j=1;
    p=P(1,:);
    p=RecursiveSearch(p,P,Psteps,j,nDims,activeDims);
    p=pv;  % get the best value.
 
else
    % P
    p=Simplex('init',P,Psteps,Psteps>0);
    %%
    % while iters<10 || ~(done || iters>niters)
    while ~(done || iters>niters)
        [mx,pt]=DoTransSearch(p,0);
        pTrans=pt(5:6);
        %%
        p=Simplex(-mx);
        done=Simplex('converged',epsi*abs(mx));
        iters=iters+1;
    end;
    p=Simplex('centroid');
    p(5:6)=pTrans;
end;
T=meMakeAffineMatrix(p);


    function [mx,p]=DoTransSearch(p,bfMode)
%         in bfMode we show only improved values.
        T=meMakeAffineMatrix(p);
        m2s=AffineTransform(m2,T);
        ccd=fftshift(real(ifftn(fc1.*conj(fftn(m2s)))));
        [mx,xs,ys]=max2di(ccd);
        p(5)=double(xs-(ndx/2+1))/ndx;
        p(6)=double(ys-(ndy/2+1))/ndy;
        better=mx>mxv;
        if better
            mxv=mx;
            pv=p;
        end;
        if displayOn && ((better&&bfMode) || (~bfMode && mod(iters,5)==0))  % show the cc function if it's better than previous
            ccalign=real(ifftn(fftn(ccd).*FourierShift(nd,-p(5:6).*nd)));
            cc=Crop(ccalign,nds);
            imacs(cc);  % show it with precise shifting
            title(num2str([iters p(5)*ndx p(6)*ndy]));
            pstr=sprintf('%8.5f   ',p(1:4));
            xlabel(pstr);
            drawnow;
        else
            cc=ccd;
        end;        
    end

    function p=RecursiveSearch(p,P,Psteps,j,jmax,activeDims)
        dim=activeDims(j);
        while p(dim)<P(2,dim)
            iters=iters+1;
            p(dim)=p(dim)+Psteps(dim);
            if j<jmax
                p=RecursiveSearch(p,P,Psteps,j+1,jmax,activeDims);
            else
                [mx,p]=DoTransSearch(p,1);
% disp([p pv mx]);
                if mx>mxv
                    mxv=mx;
                    pv=p;
                end;
            end;
        end;
        p(dim)=P(1,dim);
    end
end