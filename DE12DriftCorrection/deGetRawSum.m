function [m nim]=deGetRawSum(folderName,skipImages)
% function [sumImage nim]=deGetRawSum(folderName,skipImages)
% Read a stack of DE-12 raw images given the folder that contains them and
% compute the sum. Also write the file 'isum.mat' into the directory.

names=deGetRawImageNames(folderName);
nim=numel(names)-skipImages;
m=0;
for i=skipImages+1:numel(names)
     m1=(imread(names{i}));
     m=m+single(m1);
     imagesc(BinImage(m,4));
     title(names{i},'interpreter','none');
     drawnow;
end;
save([AddSlash(folderName) 'isum.mat'],'m','nim');
