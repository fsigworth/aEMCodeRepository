% MultiExpFit
% fit multiple exponentials to data f,x

load biaCore.mat
j=6; % 6 curves total

% dissociation
start=3636;
stop=size(biaX,1);

% association part of curve
% start=620;
% stop=3600;
% stop=2000;

q=diff(biaYs(start:stop,j));
mask=abs(q)>.05;
x=biaX(start+1:stop)-biaX(start);
q(mask)=0;
y=cumsum(q)+biaYs(start,j);

%%

alphas=[.01 .001 .0001];
nAlphas=numel(alphas);
hasConst=1;
nx=numel(x);
nf=nAlphas+hasConst;
nIters=200;

fcns=zeros(nx,nf);
newAlphas=Simplex('init',alphas);
for i=1:nIters
    alphas=newAlphas;
    for j=1:nAlphas
        fcns(:,j)=exp(-alphas(j).*x);
    end;
    if hasConst
        fcns(:,nAlphas+1)=ones(nx,1);
    end;
    a=LinLeastSquares(fcns,y);
    estY=fcns*a;
    yDiff=estY-y;
    err=yDiff'*yDiff;
    newAlphas=Simplex(err);
    if mod(i,10)==0
% %         plot(x,[y estY biaYs(start+1:stop,1:5)],'linewidth',2);
%         plot(x,[y estY biaYs(start+1:stop,1:5)],'linewidth',1.5);
        plot(x,y,'linewidth',2);
        hold on;
        plot(x,estY,'k-','linewidth',1);
        hold off;
        title(['taus: ' num2str(1./alphas,4) '     amps: ' num2str(a',4)]);
    	pause(0.1)
    end;
end;
