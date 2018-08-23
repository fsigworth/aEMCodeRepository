% Class3DConvergence.m
%  Separate out class #2 from images.

pattern='_roi.mat';
cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon96a2v');
d=dir;
bins=0:.01:1;
nb=numel(bins);
h=zeros(nb,1);
j=0;
for i=1:numel(d)
    q=strfind(d(i).name,pattern);
    disp(d(i).name);
    if numel(q)>0
        load(d(i).name);
        j=j+1;
        h(:,j)=hist(roi.pVols(1,:),bins);
        semilogy(bins,h);
        drawnow;
    end;
end;
semilogy(bins,h);
%%

goodImgs=roi.pVols(1,:)<.5;
numGoodImages=sum(goodImgs)
totalImages=numel(goodImgs)

