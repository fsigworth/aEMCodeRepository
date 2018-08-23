% re2DMakeJpegs.m
% Create graphics files in case they weren't created by re2DClassify3
nExponent=0;
frcExponent=.5;
startIter=1;
numIters=inf;

riName='ri.mat';
pa=uigetdir('.','Select a Recon directory');
if isnumeric(pa)
    return
end;
cd(pa);
[rootDir,localDir]=ParsePath(pa);
if ~exist(riName,'file')
    display('No ri.mat file found');
    return
end;

CheckAndMakeDir('jpeg/');
load(riName);
figure(3);
set(gcf,'position',[500,500,630,600]);
% Create jpegs from classes
disp(localDir);
disp('First file:');
firstIter=true;
frcTrend=zeros(ri.nCurrent/2,ri.nIters);
%%
for iter=startIter:min(numIters,ri.nIters)
    roiName=[reGetNameCode(ri,iter,1) 'roi.mat'];
            moiName=[reGetNameCode(ri,iter,1) 'moi.mat'];
    if exist(roiName,'file') && exist(moiName,'file')
        disp(roiName);
        load(roiName);
        xRoiName=[reGetNameCode(ri,iter,2) 'roi.mat'];
        if exist(xRoiName,'file') % we have both twins
            xRoi=load(xRoiName);
            xRoi=xRoi.roi;
            %        Compute the frcs
            %             disp('Computing FRCs');
            n=size(roi.classMeans,1);
            nCls=size(roi.classMeans,3);
            msk=repmat(fuzzymask(n,2,n*0.3,n*.2),1,1,nCls);
            frc=FRC(msk.*roi.classMeans,msk.*xRoi.classMeans);
            h=reFitErf(frc);
            h=max(0,(h-h(end))/(h(1)-h(end))).^frcExponent;  % normalize the filter
            hRot=ToRect(h);
            k0=1e-3;
            k1=n^nExponent;
            kx=k0+(1-hRot)./(hRot+k0);
            kx=repmat(kx,1,1,nCls)*k1;
            roi.classMeans=roi.classMeans+xRoi.classMeans;
            roi.classNorms=roi.classNorms+xRoi.classNorms;
            frcTrend(1:n/2,iter)=frc;
            figure(4);
            fs=(0:n/2-1)/(ri.pixA*n);
            plot(fs,frcTrend,'-',fs,frc,'k-',fs,h,'k--');
        else
            kx=n^nExponent;
            disp(['Using constant k' num2str(kx)]);
        end;
        
        %     =====================Do the 2D reconstruction======================
        load(moiName);
%         fCls=fftshift2(fft2(roi.classMeans));
%         nCls=size(fCls,3);
%         %         %     Wiener normalization
%         moi.refs=real(ifft2(ifftshift2((fCls)./(kx*moi.sigmaN^2+roi.classNorms))));
%         sds=sqrt(mean(mean(moi.refs.^2)));
%         k0=.001;
%         moi.refs=moi.refs./repmat(sds+k0,n,n,1);
        figure(3);
%         ImagicDisplay2(moi.refs);
        ImagicDisplay(moi.refs,1,0);
        drawnow;
        figName=[reGetNameCode(ri,iter,1) 'cls.jpg'];
                    set(gcf,'paperpositionmode','auto');
                    print(['jpeg/' figName],'-djpeg');
                    disp(figName);
        if firstIter
            firstIter=false;
            pause
        end;
    end;
end;
classes=moi.refs;
frc=frcTrend(:,max(startIter,iter-1));
figure(5);
plot(fs,[frc frc*0]);
save('jpeg/classes.mat','classes','frc','fs');
disp('done.');