function logProbWC=reComputeLogProbWC4(nt,vR,yClick,mbnOffset,ri,moi,imgIndex)
% function logProbWC=reComputeLogProbWC4(nt,vR,yClick,mbnOffset,ri,moi,imgIndex);
% Derived from logProbWC3 but assumes moi.activeTrans is nt^2 x nImgs
% identified in moi.activeT(:,imgIndex).
% 
% Compute p(c|w) x p(w|alpha,beta).  This is the product of the click error
% Gaussian (sigmaC) with the rock + bob + geometry error anisotropic
% Gaussian.  Returns an ntc x nAlphas x nIos x nGammas x nBetas set of
% frames, with origin in center.  vR is the vesicle radius - membrane
% offset.  The fields moi.sigmaC and sigmaG are used.
% All units are in pixels and degrees.
% If this is called with nt=1, returned is a vector containing only the peak
% value of the Gaussian for each alpha and beta value.

% % test code
% nt=25;
% moi.sigmaC=1.5;
% moi.sigmaG=2;
% tA=false(nt,nt);
% tA(12:15,12:15)=true;
% % tA=true(nt,nt);
% moi.activeTrans=tA(:);
% imgIndex=1;
% 
% vR=40;
% yClick=20;
% ri.angleN=[5 18 9];
% ri.angleMin=[-40 0 0];
% ri.angleStep=[20 10 20];
% ri.alphas=ri.angleMin(1)+(0:ri.angleN(1)-1)*ri.angleStep(1);
% ri.alphasI=[ri.alphas ri.alphas+180];
% ri.isos=[0 1];
% mbnOffset=10;
% tic

% Handle active flags for translation
if ~isfield(moi,'activeTrans')
    tc=true(nt^2,1);
else
    tc=moi.activeTrans(:,imgIndex);
end;
ntc=sum(tc);  % number of active translations

% nAlphas=ri.angleN(1);
% We now use the alpha list in ri.alphas instead.
dAlpha=ri.angleStep(1);
alphas=ri.alphas;
nAlphas=numel(alphas);
isos=ri.isos;
nIsos=numel(isos);

nBetas=ri.angleN(2);
dBeta=ri.angleStep(2);
betas=ri.angleMin(2)+(0:nBetas-1)*dBeta;
nGammas=ri.angleN(3);

radDeg=pi/180;
nh=(nt-1)/2;  % nt must be odd
% nhx=ceil(nh*1.5);
% ntx=2*nhx+1;  % extended size
% ctrx=nhx+1;
% Make the rotated coordinates.
[x0T,y0T]=ndgrid(-nh:nh,-nh:nh);  % Make zero at ctr,ctr
x0=x0T(tc);
y0=y0T(tc);
X=zeros(ntc,nAlphas,nIsos,'single');
Y=zeros(ntc,nAlphas,nIsos,'single');
if nt>1
    logProb=zeros(ntc,nAlphas,nIsos,1,nBetas);
    for isi=1:nIsos % loop over iso orientations
        for ia=1:nAlphas
            alphax=alphas(ia)+180*isos(isi);
            X(:,ia,isi)=cosd(alphax)*x0+sind(alphax)*y0;
            Y(:,ia,isi)=cosd(alphax)*y0-sind(alphax)*x0;
        end;
    end;
else
    logProbWC=zeros(nAlphas,nIsos,nBetas);
end;

for isi=1:nIsos
    io=isos(isi);
    vRo=vR+sign(io-.5)*mbnOffset;  % effective radius of center of particle
    for iBeta=1:nBetas;
        beta=betas(iBeta);
        % Variances of the geometric Gaussian
        varY=moi.sigmaG^2+(vRo*cosd(beta)*dBeta*radDeg).^2/12;
        varX=moi.sigmaG^2+(vRo*sind(beta)*dAlpha*radDeg).^2/12;
        ax=1/(2*varX);
        ay=1/(2*varY);
        
        %     Get the center coords of the geometric Gaussian
        % Suppose the particle is rocked in alpha by da.  The best match with the
        % reference will occur when it is unrotated by -da.  But the un-rotation will
        % make the center of the correct reference move laterallyin x and *upward* in y.
        c=1/(2*moi.sigmaC^2);        % Click Gaussian
        bx=yClick*sind(alphas);  % this ia a vector of x positions
        b1x=ax*bx/(ax+c);        % coordinate of the product of Gaussians
        by=(vRo*sind(beta)-yClick*cosd(alphas));  % also a vector
        %         by=vRo*sind(beta)-yClick;  % also a vector
        b1y=ay*by/(ay+c);       % coordinate of the product of Gaussians
        
        logPeakWC=-ax*bx.^2+ax^2*bx.^2/(ax+c)...
            -ay*by.^2+ay^2*by.^2/(ay+c);
        
        if nt>1 % want to compute the 2D log probabilities
            
            % unclear how we should normalize this.
            logNormProb=-log(2*pi*sqrt(1/(4*(ax+c)*(ay+c)))); % norm. the product
            %     logNormProb=-log(pi*sqrt(1/(ax*ay)));  % normalize the geometric Gaussian.
            %     logNormProb=-log(2*pi*sigmaC^2);  % normalize the click Gaussian.
            
            for ia=1:nAlphas  % compute the 2D gaussians
                logProb(:,ia,isi,1,iBeta)=-(ax+c)*(X(:,ia,isi)-b1x(ia)).^2 ...
                    -(ay+c)*(Y(:,ia,isi)-b1y(ia)).^2 ...
                    +logPeakWC(ia)+logNormProb;
            end;
        else  % just assign the peak value
            logProbWC(:,isi,iBeta)=logPeakWC;
            %     exp(logProbWC)
            % plot(logProbWC);
            % plot(bx,by,'.-','markersize',10)
        end;
        
    end;
end;
if nt>1
    logProbWC=repmat(logProb,1,1,1,nGammas,1);
end;


% 
% 
% toc
% 
%     figure(3);
% plot(logProbWC(:));
% %
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