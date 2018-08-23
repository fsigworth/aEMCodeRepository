function logProbWC=reComputeLogProbWC5(nt,vR,yClick,mbnOffset,ri,moi,imgIndex)
% function logProbWC=reComputeLogProbWC5(nt,vR,yClick,mbnOffset,ri,moi,imgIndex);
% Derived from logProbWC3 but handles moi.activeTrans, moi.activeAlphas and
% moi.activeRefs to prune the domain.  ntc, nai and naRefs are the number of each.
% We assume moi.activeTrans is nt^2 x nImgs
% and moi.activeAlphasI is nAI x nImgs
% identified in moi.activeT(:,imgIndex), moi.activeAlphasI(:,imgIndex)
%
% Compute p(c|w) x p(w|alpha,beta).  This is the product of the click error
% Gaussian (sigmaC) with the rock + bob + geometry error anisotropic
% Gaussian.
% Returns an ntc x nai x (naRVs) set of
% frames.  vR is the vesicle radius - membrane
% offset.  The fields moi.sigmaC and sigmaG are used.
% All units are in pixels and degrees.
% If this is called with nt=1, returned is a vector containing only the peak
% value of the Gaussian for each alpha and beta value.

% % test code
% nt=21;
% moi.sigmaC=1.5;
% moi.sigmaG=2;
% % tA=false(nt,nt);
% tA=true(nt,nt);
% moi.activeTrans=tA(:);
% imgIndex=1;
% 
% vR=40;
% yClick=35;
% ri.angleN=[5 18 9];
% ri.angleMin=[-10 0 0];
% ri.angleStep=[5 10 20];
% als=ri.angleMin(1)+(0:ri.angleN(1)-1)*ri.angleStep(1);
% ri.alphas=[als' als'+180];
% moi.activeAlphas=false(size(ri.alphas));
% moi.activeAlphas(2:4,1:2)=true;
% moi.activeAlphas=true | moi.activeAlphas; %%%
% 
% ri.isos=[0 1];
% ri.nVols=1;
% mbnOffset=5;
% tic

% Handle active flags for translation
if ~isfield(moi,'activeTrans')
    tc=true(nt^2,1);
else
    tc=moi.activeTrans(:,imgIndex);
end;
ntc=sum(tc);  % number of active translations

% Handle flags for alphas and isos
nAlphas=size(ri.alphas,1);
nIsos=size(ri.alphas,2);
if ~isfield(moi,'activeAlphas')
    acai=true(nAlphas,nIsos); % booleans
else
    acai=moi.activeAlphas(:,:,imgIndex);  % columns are alphas for given iso
end;
acai=acai(:)';
ptai=find(acai);  % pointers to active alphasI
isos=ri.isos(floor((ptai-1)/nAlphas)+1); % active inside-out values
pta=mod(ptai-1,nAlphas)+1;             % pointers to corresponding active alphas
nai=sum(acai);  % number of alpha-iso pairs

nBetas=ri.angleN(2);
nGammas=ri.angleN(3);
nVols=ri.nVols;
if ~isfield(moi,'activeRVs')
    rvc=true(nGammas,nBetas,nVols);
else
    rvc=moi.activeRVs(:,:,:,imgIndex); % nGammas x nBetas x nVols
end;
betaFlags=any(any(rvc,1),3);  % beta exists for some gamma, volume
betaPtrs=find(betaFlags);
naBetas=sum(betaFlags);  % number of active betas
gammaCounts=shiftdim(sum(rvc(:,betaFlags,:),1),1); % matrix naBetas x nVols
naRVs=sum(gammaCounts(:));

dAlpha=ri.angleStep(1);
dBeta=ri.angleStep(2);
betas=ri.angleMin(2)+(0:nBetas-1)*dBeta;

radDeg=pi/180;
nh=(nt-1)/2;  % half-width of translation; nt must be odd

% Make the rotated coordinates.
[x0T,y0T]=ndgrid(-nh:nh,-nh:nh);  % Make zero at ctr,ctr
x0=x0T(tc);  % column vector
y0=y0T(tc);
% X=zeros(ntc,nai,'single');
% Y=zeros(ntc,nai,'single');
if nt>1
    logProb=zeros(ntc,nai,naBetas,'single'); % 1 is space for gammas
    alphasI=ri.alphas(acai);
    alphasI=alphasI(:)';
    X=x0*cosd(alphasI)+y0*sind(alphasI);
    Y=y0*cosd(alphasI)-x0*sind(alphasI);
else
    logProbWC=zeros(nai,naBetas,'single');
end;

vRos=vR+sign(isos-.5)*mbnOffset;  % effective radius of center of particle
alphas=ri.alphas(pta);      % alphas without isos
alphas=alphas(:)';
c=1/(2*moi.sigmaC^2);        % Click Gaussian
t1s=ones(ntc,1);              % column vector to replicate angle values

for iBeta=1:naBetas;
    beta=betas(betaPtrs(iBeta));
    % Variances of the geometric Gaussian
    varY=moi.sigmaG^2+(vRos*cosd(beta)*dBeta*radDeg).^2/12;
    varX=moi.sigmaG^2+(vRos*sind(beta)*dAlpha*radDeg).^2/12;
    kx=1./(2*varX);
    ky=1./(2*varY);
    
    %     Get the center coords of the geometric Gaussian
    % Suppose the particle is rocked in alpha by da.  The best match with the
    % reference will occur when it is unrotated by -da.  But the un-rotation will
    % make the center of the correct reference move laterally in x and *upward* in y.
    bx=yClick*sind(alphas);  % this ia a row vector of x positions
    b1x=kx.*bx./(kx+c);        % coordinate of the product of Gaussians
    by=(vRos*sind(beta)-yClick*cosd(alphas));  % row vector of y positions
    b1y=ky.*by./(ky+c);       % coordinate of the product of Gaussians
    
    %     logPeakWC=-ax.*bx.^2+ax.^2.*bx.^2./(ax+c)...
    %         -ay.*by.^2+ay.^2.*by.^2./(ay+c);
    logPeakWC=(b1x-bx).*kx.*bx + (b1y-by).*ky.*by;
    
    if nt>1 % want to compute the 2D log probabilities
        
        % unclear how we should normalize this.
        logNormProb=-log(2*pi*sqrt(1./(4*(kx+c).*(ky+c)))); % norm. the product
        %     logNormProb=-log(pi*sqrt(1/(ax*ay)));  % normalize the geometric Gaussian.
        %     logNormProb=-log(2*pi*sigmaC^2);  % normalize the click Gaussian.
        %         Compute the log of the 2D gaussians
        logProb(:,:,iBeta)=-(t1s*(kx+c)).*(X-t1s*b1x).^2 ...
            -(t1s*(ky+c)).*(Y-t1s*b1y).^2 ...
            +t1s*(logPeakWC+logNormProb);
        
        logProbWC=zeros(ntc,nai,naRVs,'single');
        iptr=1;
        for iVol=1:nVols
            for iBeta=1:naBetas
                ng=gammaCounts(iBeta,iVol);
                if ng>0
                    logProbWC(:,:,iptr:iptr+ng-1)...
                        =repmat(logProb(:,:,iBeta),1,1,ng);
                    iptr=iptr+ng;
                end;
            end;
        end;
        
    else  % just assign the peak value, ignoring translations
        logProbWC(:,iBeta)=logPeakWC';
    end;
    
end;

%
% Test Code
% toc
% %%
%     figure(3);
% plot(logProbWC(:));
% q=reshape(logProbWC,nt,nt,10,9,18);
% plot(squeeze(q(:,5,:,2,10)))
% legend(num2str(alphasI'))
% 
% imovie(q(:,:,:,1,10),.5)
% %
% return
% 
% 
% % test calc
% 
% %%
% varX=10;
% sigmaC=1;
% xOrg=80;
% b=xOrg;
% x=(-100:100)';
% 
% lr=-(x-b).^2/(2*varX)-x.^2/(2*sigmaC^2);
% 
% 
% 
% a=1/(2*varX);
% b=xOrg;
% c=1/(2*sigmaC^2);
% d=a+c;
% b1=a*b/d;
% lq=-d*(x-b1).^2-a*b^2+a^2*b^2/d;
% plot([lr lq]);
% 
