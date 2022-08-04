% CropScreenshots.m
outPath='../Imani2/';
doWrite=1;
figure(1);
clf
% set(gcf,'toolbar','none');
% set(gcf,'menubar','none');
fs=14;
mysubplot(1,1,1);
d=dir;
nd=numel(d);
inds=37-[3 6 12 32]; % start is 34
nSets=numel(inds)-1;
ul=[48 222
    121   105
    14   161
    1   1];
lr=[2318 1497
    2296 1335
    2436  1525
    2880 1800];
border=50;
k=0;
for iSet=1:3
    for i=inds(iSet):-1:inds(iSet+1)+1
        k=k+1;
        inName=d(i).name;
        m=imread(inName);
        mc=m(ul(iSet,2):lr(iSet,2),ul(iSet,1):lr(iSet,1),:);
        image(mc);
        text(border,lr(2)-ul(2)-border,num2str(k),'fontsize',fs);
        axis equal off;
        set(gcf,'name',inName);
        drawnow;
        outName=[outPath 'c' num2str(k) '.png'];
        disp(outName);
        % outName=[imgBaseName num2str(k) '.png'];
        if doWrite
            imwrite(mc,outName,'png');
        end;
        disp([iSet i]);
    end;
end;