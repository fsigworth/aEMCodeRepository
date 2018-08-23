function logProbWC=reComputeLogProbWC(n,sigmaC,sigmaG,vR,yClick,alphas,dAlpha,beta,dBeta)
% function logProbWC=reComputeLogProbWC(n,sigmaC,sigmaG,vR,yClick,alphas,dAlpha,beta,dBeta)
% Compute the product of the click error Gaussian (sigmaC) with the
% rock + bob + geometry error anisotropic Gaussian.  Returns an n x n x nAlphas set of
% frames, with origin in center.  vR is the vesicle radius - membrane
% offset.
% All units are in pixels and degrees.
% If this is called with n=1, returned is a vector containing only the peak
% value of the Gaussian for each alpha value.
% n=64;
% sigmaC=4;
% sigmaG=10;
% vR=40;
% yClick=35;
% alphas=(0:10:90)';
% beta=60;
% dAlpha=.1;
% dBeta=.1;
% tic


nAlphas=numel(alphas);
radDeg=pi/180;

% Variances of the geometric Gaussian
varY=sigmaG^2+(vR*cosd(beta)*dBeta*radDeg).^2/12;
varX=sigmaG^2+(vR*sind(beta)*dAlpha*radDeg).^2/12;
ax=1/(2*varX);
ay=1/(2*varY);

%     Get the center coords of the geometric Gaussian
% Suppose the particle is rocked in alpha by da.  The best match with the
% reference will occur when it is unrotated by -da.  But the un-rotation will
% make the center of the correct reference move laterallyin x and *upward* in y.
c=1/(2*sigmaC^2);        % Click Gaussian
bx=yClick*sind(alphas);  % this ia a vector
b1x=ax*bx/(ax+c);        % coordinate of the product of Gaussians
by=vR*sind(beta)-yClick*cosd(alphas);  % also a vector
b1y=ay*by/(ay+c);       % coordinate of the product of Gaussians

logPeakWC=-ax*bx.^2+ax^2*bx.^2/(ax+c)...
          -ay*by.^2+ay^2*by.^2/(ay+c);
if n>1 % want to compute the 2D log probabilities
    ctr=floor(n/2+1);
    [x,y]=ndgrid(1-ctr:n-ctr,1-ctr:n-ctr);  % Make zero at ctr,ctr

    % unclear how we should normalize this.
     logNormProb=-log(2*pi*sqrt(1/(4*(ax+c)*(ay+c)))); % norm. the product
%     logNormProb=-log(pi*sqrt(1/(ax*ay)));  % normalize the geometric Gaussian.
%     logNormProb=-log(2*pi*sigmaC^2);  % normalize the click Gaussian.

    logProbWC=zeros(n,n,nAlphas);
    for i=1:nAlphas  % compute the 2D gaussians
        logProbWC(:,:,i)=-(ax+c)*(x-b1x(i)).^2-(ay+c)*(y-b1y(i)).^2 ...
                         +logPeakWC(i)+logNormProb;
    end;
% q=sum(sum(exp(logProbWC)));
% % semilogy(alphas,q(:))
% q(:)
%     q=max(logProbWC(:));
%     for i=1:nAlphas
% %         imac(256+5*logProbWC(:,:,i));
%         imac(256*exp(logProbWC(:,:,i)-q));
%         pause(0.1);
%     end;
else  % just assign the peak value
    logProbWC=logPeakWC;
%     exp(logProbWC)
% plot(logProbWC);
% plot(bx,by,'.-','markersize',10)
end;
% toc

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