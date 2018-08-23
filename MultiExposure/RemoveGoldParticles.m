% function mOut=RemoveGoldParticles(m,fcs)

% rootPath='/Volumes/TetraData/EMWork/Hideki/DNA-LND/data/FigA/Micrograph/';
% d=dir(rootPath);
% nd=numel(d);
% k=3;
% names={};
% for i=1:nd-k
%     names{i}=d(i+k).name;
% end;
%
% pixA=3.9;
% m=meReadImages(names,16, 0, pixA);
%
%%
[nx ny nim]=size(m);
if nargin<3
    fcs=[.04 .02 .01];
end
fcs(end+1:nim)=.01;  % default extension
nsd=1;
niter=4;

nim=1;

mOut=m;
for i=1:nim
    m1=m(:,:,i);
    m1=m1-mean(m1(:));
    m1f=GaussFilt(m1,.2);
    s=std(m1f(:));
    msk=m1f<-nsd*s;
    mskx=GaussFilt(msk,fcs(i))<.01;
    m2=m1.*mskx;
    for j=1:niter
        m2f=GaussFilt(m2,fcs(i)/10);
        m2=m1.*mskx+m2f.*(1-mskx);
        imacs(Crop(m2,4096));
        title(j);
        drawnow; 
    end;
    mOut(:,:,i)=m2;
    
end;



