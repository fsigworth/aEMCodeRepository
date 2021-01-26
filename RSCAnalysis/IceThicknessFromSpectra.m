% IceThicknessFromSpectra.m


movieDir='Movies/';
pixA=1.13;
spBin=1;
mvNames=cell(0,1);
nmv=0;
d=dir(movieDir);
for i=1:numel(d)
    [pa,nm,ex]=fileparts(d(i).name);
    if ~d(i).isdir && strcmp(ex,'.tif')
        nmv=nmv+1;
        mvNames{nmv,1}=[movieDir d(i).name];
    end;
end;

for i=1:nmv
    disp(mvNames{i});
    mv=single(ReadMovie(mvNames{i}));
    %
    if i==1
        n=size(mv);
        nFrames=n(3);
        n(3)=[];
        nPad=NextNiceNumber(n,5,8)
        nsp=min(nPad)/(2*spBin);
        sp1s=zeros(nsp,nmv,'single');
        df=spBin/(nsp*spBin*2*pixA);
        freqs=(1:nsp)*df;
    end;
    %
    for k=1:nFrames
        mvf=mv(:,:,k);
        mvfc=Crop(mvf,nPad,0,median(mvf(:)));
        sp1s(:,i)=sp1s(:,i)+RadialPowerSpectrum(mvfc,0,spBin);
        fprintf('.');
    end;
    sp1s(1:7,i)=sp1s(8,i);
    fprintf('\n');
    plot(freqs,sp1s(:,1:i));
    drawnow;
end;
