function xytfr=MRAlign(imgs,refs,mode,thetas,doFlip,transSD,debug)
% function xytfr=MRAlign(imgs,refs,mode,thetas,doFlip,transSD,debug);
% Multi-reference alignment of an image stack imgs to a stack of references
% refs.  It is up to the user to pre-mask the refs and images, so that they
% are zero at the periphery.
% The returned xytfr matrix is nimages x 6 and gives the xshift, yshift,
% theta, flip about x (0 or 1) and reference number.  The 6th element is the maximum
% cc value.  theta is in degrees.  The modes are
% 1: translate only
% 2: rotate only
% 3: simultaneous translate + rotate
% thetas is a list of theta values to search, in degrees.
% transSD, if it's >0, is actually the ratiosigma_trans / sigma_noise,
% where sigma_trans is in pixels and sigma_noise is the pixel noise.
% defaults are: mode=3, thetas=0:10:359, doFlip=0, transSD=0, debug=0.

[n, ny, nim]=size(imgs);
[nxr, nyr, nrefs]=size(refs);

if nxr~=n  % change the refs size to match the images
    refs=Downsample(refs,n,1);
end;

if nargin<7
    debug=0;
end;
if nargin<6
    transSD=0;
end;
if transSD==0
    logPrior=0;
else
    % actual prior should be s_n^2 * |r-r0|^2/(2*s_r^2), where s_n is noise sigma,
    %     and s_r is translation sigma.  We assume s_n=1.
    logPrior=-Radius(n).^2/(2*transSD^2);
end;
if nargin < 5
    doFlip=0;
end;
nflip=1+(doFlip>0);

if nargin<4
    thetas=0:10:359;
end;
if nargin<3
    mode=3;
end;

if mode==1  % translate only
    nthetas=1;
    thetas=0;
else
    nthetas=numel(thetas);
    thetasW=[thetas(end)-360; thetas(:); thetas(1)+360]; % wrap around
end;

if debug
    nrefs
    nthetas
    disp('Rotating references...');
end;
reflib=single(zeros(n,n,nthetas,nflip,nrefs));

for j=1:nrefs
    ref=refs(:,:,j);
    ref=ref/sqrt(ref(:)'*ref(:));  % normalize each reference
    reflib(:,:,:,1,j)=rsRotateImage(ref,thetas); % make set of rotations
    if nflip>1
        for k=1:nthetas
            % flipud flips the first (x) coordinate, but doesn't move the fftcenter to
            % the right place if nx is even.  We flip the rotated reference.
            flipref=circshift(flipud(reflib(:,:,k,1,j)),[mod(n+1,2) 0]);
            reflib(:,:,k,2,j)=flipref;
        end;
    end;
end;

if mode ~=2  % We need fts for trans alignment
    freflib=conj(fft2(reflib));
end;

if debug
    disp('Aligning...');
end;

% Do the alignment
% xytfr=zeros(nim,5);

ctr=ceil((n+1)/2);  % image center coordinate

xytfr=zeros(nim,6);

% cc=zeros(nthetas,1);
% ccl=zeros(n,n,nthetas+2);
% sx=zeros(nflip,1);
% sy=zeros(nflip,1);
% sj=ones(nflip,1);
% mx=zeros(nflip,1);
% mxl=zeros(nrefs,1);
% tx=zeros(nrefs,1);
% ty=zeros(nrefs,1);
% tj=zeros(nrefs,1);
% tfl=zeros(nrefs,1);
%     Compute the inner product of each image with each reference
% parfor i=1:nim
for i=1:nim
    cc=zeros(nthetas,1);
    ccl=zeros(n,n,nthetas+2);
    sx=zeros(nflip,1);
    sy=zeros(nflip,1);
    sj=ones(nflip,1);
    mx=zeros(nflip,1);
    mxl=zeros(nrefs,1);
    tx=zeros(nrefs,1);
    ty=zeros(nrefs,1);
    tj=zeros(nrefs,1);
    tfl=zeros(nrefs,1);
    im=imgs(:,:,i);
    if mode ~=2
        fim=fft2(im);
    end;
    for j=1:nrefs
        for l=1:nflip
            if mode==2
                for k=1:nthetas
                    rf=reflib(:,:,k,l,j);
                    cc(k)=im(:)'*rf(:);  % just the inner product
                    [mx(l), sj(l)]=max1di(cc);  % best theta for each ref
                end;
            else
                for k=1:nthetas % compute 3D cross-correlation
                    ccl(:,:,k+1)=fftshift(real(ifftn(freflib(:,:,k,l,j).*fim)))+logPrior;
                end;
                if nthetas>1
                    ccl(:,:,1)=ccl(:,:,nthetas+1);  % wrap lower end
                    ccl(:,:,nthetas+2)=ccl(:,:,1);  % wrap upper end
                    [mx(l), coords]=max3di(ccl);
                    sx(l)=coords(1);
                    sy(l)=coords(2);
                    sj(l)=max(1,min(nthetas+1,coords(3))); % force in bounds
                else
                    [mx(l), sx(l), sy(l)]=max2di(ccl(:,:,2));
                end;
            end;
        end;  % for l
        [mxl(j), ifl]=max(mx);  % Get the best over mirror images
        tx(j)=sx(ifl);
        ty(j)=sy(ifl);
        tj(j)=sj(ifl);
        tfl(j)=ifl;
    end;
    [mxcc, iref]=max(mxl);  % best ref overall
    uj=tj(iref);
    theta=interp1(thetasW,uj); % thetas1 wraps around
    xytfr(i,:)=[tx(iref)-ctr ty(iref)-ctr theta tfl(iref)-1 iref mxcc];
    
%     if mod(i,1000)==0
%         disp(i);
%     end;
end;  % loop over images i

% else
%
%     % shift or shift-and-rotate.  The former simply has nthetas=1.
%     ct=n/2+1;
%     if debug
%         disp('Rot and shift determination...');
%     end;
%     for i=1:nim
%         cc1=single(zeros(n,n,nthetas));  % translational cc
%         iflip=zeros(nrefs,1);
%         coords=zeros(3,nrefs);
%         mx=zeros(nflip,1);
%         cmx=zeros(3,nflip);
%         mxvals=zeros(nrefs,1);
%         im=imgs(:,:,i);
%         fim=fftn(im);
%         for j=1:nrefs
%             for l=1:nflip
%                 for k=1:nthetas
%                     cc1(:,:,k)=fftshift(real(ifftn(freflib(:,:,k,l,j).*fim)))+prior;
%                 end;
%                 if nthetas>1
%                     [mx(l), cmx(:,l)]=max3di(cc1);
%                 else
%                     [mx(l), cmx(1,l), cmx(2,l)]=max2di(cc1);
%                 end;
%             end;
%             [mxvals(j), iflip(j)]=max(mx);
%             coords(:,j)=cmx(:,iflip(j));  % encoded x,y,t for ref. j
%         end; % loop over references j
%         [mxcc, refj]=max(mxvals);
%         theta=dt*(coords(3,refj)-(nthetas+1)/2);
%         % q=[coords(1:2,refj)'-ct theta iflip-1 refj mxcc]
%         xytfr(i,:)=[ct-coords(1:2,refj)' theta iflip(refj)-1 refj mxcc];
%
%         if mod(i,1000)==0
%             disp(i);
%         end;
%     end;  % loop over images i
%     end;
%
%
% % if nargout>1  % we return aligned images
% if debug
%     disp('Rotating output images...');
% end;
% alignedimages=zeros(n,n,nim);
% for j=1:nim
%     im=imgs(:,:,j);
%     theta=-xytfr(j,3);
%     iflip=xytfr(j,4);
%     if intShifts
%         im=circshift(imgs(:,:,j),-round(xytfr(j,1:2)));
%     else
%         im=shiftf(imgs(:,:,j),-xytfr(j,1:2));
%     end;
%     if theta ~= 0
%         im=grotate(im,theta);
%     end;
%     if iflip
%         im=circshift(flipud(im),[mod(n+1,2) 0]);  % flip after rotate and shift.
%     end;
%     alignedimages(:,:,j)=im;
% end;
%
