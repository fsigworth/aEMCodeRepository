%  SumAllImages.m

[names,pa]=uiGetFilenames;
cd(pa);

%%
for i=1:numel(names)
    name=names{i};
    disp(name);
%     m=ReadEMFile(name);
    m=sum(single(ReadMovie(name)),3);
    if i==1
        msum=m;
    else
        msum=msum+m;
    end;
end;
%%

msumo=RemoveOutliers(msum);
%%
figure(1);
SetGrayscale;
subplot(1,2,1);
imacs(GaussFilt(msum,.1));

% %%
[refName,pa]=uiGetFilenames;
mr=ReadEMFile([pa refName{1}]);
subplot(1,2,2)
mr=RemoveOutliers(mr);
%

% mrc=Crop(BinImage(mr,2),1920,0,mean(mr(:)));
mrc=Crop(BinImage(mr,2),3840,0,mean(mr(:)));
% imacs(GaussFilt(msum./max(mrc,.9),.1));
% imacs(GaussFilt(mrc.*msum,.1));
imacs(GaussFilt(mrc,.1));
