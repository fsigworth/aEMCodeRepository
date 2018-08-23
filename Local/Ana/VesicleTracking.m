% VesicleTracking.m

cd('/Users/fred/Box Sync/VesicleMovies')

load modelVesicle.mat
load Vesicle23.mat

nRep=100;
n=256;
ves=Crop(Downsample(modelVesicle,384),n)*8;
%   scale up by 8 allows good subtraction
mvd=Downsample(mvx,384,1);
mvd=Crop(mvd,n,1);
nf=17;  % number of frames to use
mvd=mvd(:,:,1:nf);  % just the low-defocus part.

f0Vals=[.03 .04 .05 .06 .07];
filtExponent=0;
sigmaRef=1;

nf0=numel(f0Vals);
varMeans=zeros(nf0,1);

ref=ves+sigmaRef*randn(n,n);

for if0=1:nf0
    f0=f0Vals(if0);
    
    freqs=RadiusNorm(n);
    % H=1./((freqs/f0).^2+1);
    H=abs(freqs).^filtExponent./((freqs/f0).^2+1);
    
    varVals=zeros(nRep,1);
    for j=1:nRep
        
        % fake data
        sigma=15;
        mvd=sigma*randn(n,n,nf)+repmat(ves,1,1,nf);
        
        mx=sum(mvd,3);
        
        cc=zeros(n,n,nf);
        xs=zeros(nf,1);
        ys=zeros(nf,1);
        for k=1:nf
            cc(:,:,k)=fftshift(ifftn(...
                ifftshift(H).*fftn(mvd(:,:,k)).*conj(fftn(ref))  ));
            [mxv,xs(k),ys(k)]=max2di(cc(:,:,k));
        end;
        xs=xs-(n/2+1);
        ys=ys-(n/2+1);
        
        var=(xs'*xs+ys'*ys)/nf;
        
        varVals(j)=var;
        
        subplot(2,3,1);
        plot(xs,ys);
        title(sqrt(var));
        subplot(2,3,2);
        plot([xs ys])
        title(f0);
        subplot(2,3,3);
        imags(mx);
        title(j);
        subplot(2,3,4);
        imags(H);
        subplot(2,3,5);
        imags(Crop(cc(:,:,1),32));
        drawnow;
    end;
    varMeans(if0)=median(varVals);
    subplot(2,3,3);
    hist(varVals);
end;
%%
subplot(2,3,6);
plot(varMeans);
ylabel('rms error');
xTickLabel=cell(0,0);
for if0=1:nf0
    xTickLabel{if0}=num2str(f0Vals(if0));
end;
set(gca,'xTickLabel',xTickLabel);
xlabel('f0');
title(['filtExponent = ' num2str(filtExponent)]);
return

%%

imovie(Crop(cc,64,1),.1);
