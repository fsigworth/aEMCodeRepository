function [mergeMats, corrCTPars, mergeSNR]=meAlignMultiExposures(m,pixA,rawCTFs,doses,ds,oldMats,showGraphics)
% function [mergeMats corrCTPars]=meAlignMultiExposures(m,pixA,CTFitPars,doses,ds,oldMats,showGraphics)
% Given the stack of images m (n x n x nim), align them all, returning the
% transformation matrices mergeMats (3 x 3 x nim).  Alignments are carried
% out sequentially (m(:,:,2) aligned to m(:,:,1), m(:,:,3) aligned to
% m(:,:,2) etc.) and then the products of the transformations are taken.
% In the end, if you compute
% ma(:,:,j)=AffineTransform(m(:,:,j), mergeMats(:,:,j)) for each j, the
% resulting ma images will all be aligned.
% pixA is the angstroms/pixel.  CTFitPars is an array (nim) of the CTF
% parameter structures. ds is the optional downsample factor for aligning;
% default is 4.
% A fixed noise correction is performed.  A cross-correlation is computed
% between the first and last image.  It is assumed that the defocus
% difference between these will be large enough that the actual CC peak
% will be separate from the zero-lag peak.  We take the zero-lag peak to be
% spurious (from common gain references) and we null that peak out by
% adding appropriate anti-correlated noise.
%
% If oldMats is given and has sufficient elements, the CTPars are corrected
% and oldMats is returned as mergeMats.
%
% Statistics from 14000 mi.mergeMatrices show that T(1,1)--roughly the mag
% is between
% 1 and 1.08; T(1,2)--roughly theta-- is between -.01 and .002. Thus a brute-force search in
% these ranges is justified.  fs 25 May 14
%
% Variable naming conventions:
%   d means downsampled
%   z means the zero-frequency point is in the center.
antiNoiseFactor=1;  % should be 1, I don't know why not..

mergeSNR=1;
if nargin<5
    ds=4;  % downsampling factor for alignment, relative to original micrograph
end;
ds1=4;  % further shrink by 2

if nargin<6
    oldMats=[];
end;

if nargin<7
    showGraphics=0;
end;

% if nargin>5 && any(Ps(:))
%     refineOnly=1;
% else
refineOnly=0;
%     Ps=[];
disp('Aligning exposures for merging');
% end;

doFixedNoiseCorrection=1;
displayFNC=0;  % Don't display it.
windowFractionalWidth=16;  % used to be 128, but this seems better.


%%
[nx, ny, nim]=size(m);

if numel(oldMats)<9*nim % not a complete set of matrices, need to fit
    
    n0=[nx ny];
    nxd=nx/ds;
    nyd=ny/ds;
    nd=[nxd nyd];
    
    pixAd=pixA*ds;  % Angstroms per pixel
    
    win=single(SquareWindow(n0,nx/windowFractionalWidth));
    winFraction=sum(win(:))/numel(win);
    % fy=fftn(y.*win);  % FT of compensation noise
    
    
    nit=[0 40 200];  % iterations of the Simplex fitting.
    
    
    % allocate arrays
    
    Tmats=zeros(3,3,nim);     % transform matrices
    Tmats(:,:,1)=eye(3);
    md=single(zeros(nxd,nyd,nim));  % downsampled images
    ctz=single(zeros(nxd,nyd,nim)); % corresponding CTFs
    fd=complex(single(zeros(nxd,nyd,nim)));  % FTs
    fd0=fd;
    dsmaskz=single(zeros(nx,ny,2));  % downsampling masks
    dsmaskz(:,:,1)=single(fuzzymask(n0,2,nd*.45,nd*.05));
    dsmaskz(:,:,2)=single(fuzzymask(n0,2,nd*.225,nd*.05));
    
    % Do fixed-noise correction on the two highest-defocus images.
    if nim>1 && doFixedNoiseCorrection
        [fy, cc]=meMakeCompNoiseAuto(m(:,:,nim-1),m(:,:,nim));
        dose0=sqrt(prod(doses(nim-1:nim)));
    else
        fy=0;
        dose0=1;
    end;
    
    % Process and downsample the images.
    ndis=64;
    %%
    for i=1:nim
        f0=fftn(m(:,:,i).*win);
        f1=f0+(-1)^i*doses(i)/dose0*fy*winFraction*antiNoiseFactor;  % put in fixed-noise correction
        %     make the FT of downsampled images with anti-correlation constant
        fd(:,:,i) =ifftshift(Crop(fftshift(f1).*dsmaskz(:,:,1+(i>1)),nd));
        
        % Optional: display the residual correlations to check fixed-noise.
        if showGraphics && doFixedNoiseCorrection && displayFNC
            figure(1); clf;
            SetGrayscale;
            %         same as fd, but with no anti-noise
            fd0(:,:,i)=ifftshift(Crop(fftshift(f0).*dsmaskz(:,:,1+(i>1)),nd));
            if i>1
                ccs0=Crop(fftshift(real(ifftn(fd0(:,:,i-1).*conj(fd0(:,:,i))))),ndis);
                ccs=Crop(fftshift(real(ifftn(fd(:,:,i-1).*conj(fd(:,:,i))))),ndis);
                subplot(3,nim,i+nim);  % second row of figure, column i
                imacs(ccs0);
                title('CC not corrected');
                drawnow;
                subplot(3,nim,i+2*nim); % 3rd row of figure
                imacs(ccs);
                title('CC corrected');
                drawnow;
            end;
        end;
        % Construct the downsampled images
        md(:,:,i)=real(ifftn(fd(:,:,i)));
        if showGraphics && displayFNC
            subplot(3,nim,i);  % top row of figure
            imacs(md(:,:,i));  % original image
            axis off
            title(['Image ' num2str(i)]);
            drawnow;
        end;
        ctz(:,:,i)=single(CTF(nd,pixAd,rawCTFs(i)));
    end;
    fdz=RadiusNorm(nd);
    %%
    for i=2:nim  % find Tmat such that Tmat(image i) matches (image i-1)
        % Matching filter.  Empirically a sqrt(c1*c2) filter looks good.
        %     filt=ifftshift( sign(ctz(:,:,i-1)).*sign(ctz(:,:,i))...
        %         .*sqrt( abs(ctz(:,:,i-1).*ctz(:,:,i)).*fdz ) );  % added sqrt(fdz)
        
        filt=sign(ctz(:,:,i-1)).*sign(ctz(:,:,i))...
            .*( abs(ctz(:,:,i-1).*ctz(:,:,i)).*fdz );  % added fdz factor
        
        %     disp(' aligning');
        %     disp(' displayed parameters are: theta (radians), mag, stx, sty')
        if ~refineOnly  % Make a de novo fit.
            
            %     First, do a brute-force alignments at low resolution
            md1=Downsample(md,nd/ds1,1);
            filt1=Crop(filt,nd/ds1);
            %          for brute force
            PB=[-.03 .8 0 0 0 0;       % +/- 2 degrees, 20% mag changes
                .03 1.2 0 0 0 0];
            PstepsB=[.0071 .0071 0 0 0 0];  % 1/140 radian, mag change.
            %
            [P,T,cc01]=meMergeAligner3(md1(:,:,i-1),md1(:,:,i),filt1,PB,PstepsB,0,showGraphics);
            
            mx1=max(cc01(:))/std(cc01(:));
            mergeSNR(i)=mx1;
            disp(['Initial alignment SNR = ' num2str(mx1)]);
            
            %     %% Optimize the alignment
            % First, vary stretch only
            P(5:6)=0;  % Get rid of translation
            Psteps2=[0 0 10/nx 10/ny 0 0];
            P=meMergeAligner3(md(:,:,i-1),md(:,:,i),filt,P,Psteps2,nit(2),showGraphics);
            
        end;
        %         P=Ps(i,:);
        Psteps1=[.002 .002   0    0    0 0];
        Psteps2=[  0   0  .0002 .0002  0 0];
        
        % Finally, vary all parameters
        P(5:6)=0;
        Psteps3=(Psteps1+Psteps2);
        [P, T]=meMergeAligner3(md(:,:,i-1),md(:,:,i),filt,P,Psteps3,nit(3),showGraphics);
        % Store the transformation matrix.
        Tmats(:,:,i)=T;
    end;
    
    % Form the composite matrices.  to match image 3 to image 1, we operate
    % with T1*T2*T3, while image 1 is just operated on by T1.
    mergeMats=Tmats;  % Allocate it, and copy the first matrix.

for i=2:nim % form the composite matrices
    mergeMats(:,:,i)=mergeMats(:,:,i-1)*Tmats(:,:,i);
end;

else
    mergeMats=oldMats;
    mergeSNR=0;
end;  % if numel(oldMats) sufficient.

% Update the CTF parameters
% in view of the magnification changes.
corrCTPars=rawCTFs;  % Most fields are unchanged.
for i=2:nim
    q=mergeMats(1:2,1:2,i);   % Get the rotate/mag part
    scalesq=prod(sqrt(sum(q.^2)));  % magnification squared
    corrCTPars(i).defocus=rawCTFs(i).defocus*scalesq;
    corrCTPars(i).deltadef=rawCTFs(i).deltadef*scalesq;
end;

