function deComparePowerSpectra(sp1,sp2,p)
% display zero-centered spectra

fsize=12; 

[nx0 ny0]=size(sp1);
    n0=[nx0 ny0];
     p.imageCropFraction=1;
%     p.spectCropFraction=.7;
%     p.upscale1=3;
%     p.upscale2=3;
%     p.displayExp2=.2;
%     p.displayExp1=.3;
%     p.kStart=-3.5;
%     p.fexp1=0.4;
%     p.displayFc=.01;
%     p.lfbin=100;  % number of lf spectral points to ignore.
        spectCropSize=round([nx0 ny0]*p.imageCropFraction*p.spectCropFraction);
    f=RadiusNorm(spectCropSize)*p.spectCropFraction;

    startSp=Crop(sp1,spectCropSize).*exp(p.k1*f.^p.fexp1);
    [ndx ndy]=size(startSp);
    fsx=(-ndx/2:ndx/2-1)/ndx*p.spectCropFraction;
    fsy=(-ndy/2:ndy/2-1)/ndy*p.spectCropFraction;
    spSign=sign(startSp);
    startSp=abs(startSp);
%     startSp=max(0,startSp);
    dissp=GaussFilt(spSign.*startSp.^p.displayExp1,p.displayFc);
    startSpr=Radial2((sp1));
    mxv=startSpr(round(numel(startSpr)*p.satFreq1)+1);
    mnv=startSpr(round(numel(startSpr)*.4)+1);
    n1DSpect=spectCropSize(2)/2;
    mnv=mean(startSpr(round(n1DSpect/2)+1:n1DSpect));
    mscale=256/((mxv-mnv));
    subplot(2,2,1);
    imac(fsx,fsy,(dissp-mnv)*mscale+p.offs1);
    set(gca,'fontsize',fsize);
    title('A. Original');
    xlabel('Spatial frequency, pix^{-1}');
    drawnow;
    
    endSp=Crop(sp2,spectCropSize);
%     endSp=max(endSp,0);
%     diesp=GaussFilt(endSp.^p.displayExp2,p.displayFc);
    spSign=sign(endSp);
    endSp=abs(endSp);
%     startSp=max(0,startSp);
    diesp=GaussFilt(spSign.*endSp.^p.displayExp2,p.displayFc);
    endSpr=Radial2((sp2));

        mxv=endSpr(round(numel(endSpr)*p.satFreq2)+1);
    mnv=mean(endSpr(round(n1DSpect/2)+1:n1DSpect));
    mscale=256/((mxv-mnv));
    subplot(2,2,2);
    imac(fsx,fsy,(diesp-mnv)*mscale+p.offs2);

    
    set(gca,'fontsize',fsize);
    title('B. Compensated');
    xlabel('Spatial frequency, pix^{-1}');

    drawnow;
    disp('Computing radial spectra');
    
    subplot(223);
    startSpr=startSpr/prod(n0);
    n1DSpect=spectCropSize(2)*0.5;
    fs=(1:n1DSpect)/n1DSpect/2*p.spectCropFraction;
    startSpr=startSpr(1:n1DSpect);
%     startSpr=startSpr.*exp(p.k1*fs.^p.fexp1)';
    startSpr(1:p.lfbin)=startSpr(p.lfbin);
    plot(fs,startSpr(1:n1DSpect));
    smin=min(startSpr);
    smx=startSpr(p.lfbin);
    axis([0 inf smin-0.1*(smx-smin) smx]);
    xlabel('Spatial frequency, pix^{-1}');
    ylabel('Spectral density');
    set(gca,'fontsize',fsize);
    title('C');
    xlabel('Spatial frequency, pix^{-1}');
    
    subplot(224);
    endSpr=endSpr/prod(n0);
    endSpr(1:p.lfbin)=endSpr(p.lfbin);
    plot(fs,endSpr(1:n1DSpect));
    smin=min(endSpr(1:n1DSpect));
    smx=endSpr(p.lfbin);
    axis([0 inf smin-0.1*(smx-smin) smx]);
    set(gca,'fontsize',fsize);
    title('D');
    xlabel('Spatial frequency, pix^{-1}');
    ylabel('Spectral density');
    drawnow;
    
    
