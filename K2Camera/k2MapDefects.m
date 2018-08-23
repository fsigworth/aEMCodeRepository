% MapK2Defects.m
% cd('/Users/fred/EMWork/Hideki/151216_0/movie-frames/1');

load('defectsgrid50001.mat');
df=defectsgrid50001(:);
df(isnan(df))=[];
dfx=floor(df(1:2:end)/2)+1;
dfy=floor(df(2:2:end)/2)+1;
np=numel(dfx);
map=zeros(3840,3840,'single');
for i=1:np
    map(dfx(i),dfy(i))=map(dfx(i),dfy(i))+1;
end;
%%
imags(GaussFilt(map,.02));

%%
m=ReadMovie('grid5_Dec16_16.11.26.tif');
ms=sum(single(m),3);
%%
vals=1:1000;
h=hist(ms(:),1:1000);
semilogy(vals,h);
%%
thresh=median(ms(:))+6*std(ms(:));
map=ms>thresh;
map=sqrt(max(0,ms-thresh));
imags(GaussFilt(map,.02));
title(sum(map(:)));

%%
[mx,mapd]=RemoveOutliers(ms);
sum(mapd(:))
imags(GaussFilt(mapd,.02));
