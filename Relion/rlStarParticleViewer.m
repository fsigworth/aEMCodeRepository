% rlStarParticleViewer.m
% Starting in the experiment directory, open a data.star file from a 2D or 
% 3D classification, and step through micrographs, marking particles on them.

disp('Getting a data.star file');
[sName,sPath]=uigetfile('*.star');
if isnumeric(sPath) % user clicked Cancel
    return
end;
starName=[sPath sName];
[names,dat]=ReadStarFile(starName);

d=dat{1};
nLines=numel(d.rlnCoordinateX);
disp([num2str(nLines) ' particles.']);

%%
p=struct;
p.ds=4;
p.boxRad=25;
p.nst=8;
jpegDir='jpeg/'; % place to put jpegs if we are writing them.

[micNames,p.iMics,p.iParts]=unique(d.rlnMicrographName);
nMics=numel(micNames);
nClasses=max(d.rlnClassNumber);
% clsBools=cell(nClasses,1);
clsBools=cell(nClasses,1);
for i=1:nClasses
    clsBools{i}=d.rlnClassNumber==i;
end;

% Make the defocus histogram
defs=(d.rlnDefocusU+d.rlnDefocusV)/2e4;
p.micDefs=defs(iMics);
xs=min(defs):.1:max(defs);
h=zeros(numel(xs),nClasses);
for i=1:nClasses
    h(:,i)=hist(defs(clsBools{i}),xs)/sum(clsBools{i});
end;
figure(1);
bar(xs,h);
xlabel('Defocus, \mum');
ylabel('Relative frequency');
    
%% Mark particles on micrographs
% Box colors for the classes
p.colors=[0 1 0
       .8 .3 0
       .8 .2  .3];

figure(2);

ind=1;
ShowData(ind,d,p);
b=' ';
while b~='q'

    switch b
        case 'n'
            if ind>=nMics
                beep
            else
                ind=ind+1;
                ShowData(ind,d,p);
            end;
        case 'p'
            if ind<2
                beep
            else
                ind=ind-1;
                ShowData(ind,d,p);
            end;
        case 'R' % Run to the end, saving each image as a jpeg.
%                  This can't be interrupted by a keystroke.
            name=d.rlnMicrographName{p.iMics(ind)};
            [pa,nm,ex]=fileparts(name);
            disp(name);
            CheckAndMakeDir(jpegDir);
            print('-djpeg',[jpegDir nm '.jpg']);
            if i<nMics
                ind=ind+1;
                ShowData(ind,d,p);
            else
                break;
            end;

    end;
    if b~='R'
[px,py,b]=ginput(1);
    end;
    % b=char(b);
end;
disp('Done.');

%%
function ShowData(micIndex,d,p)
name=d.rlnMicrographName{p.iMics(micIndex)};
m=ReadMRC(name);
n0=size(m);
n=NextNiceNumber(size(m),5,8);
mc=Crop(m,n,1,mean(m(:)));
md=Downsample(mc,n/p.ds);
mf=GaussFilt(md,.1);
msc=256/(p.nst*std(mf(:)))*(mf-mean(mf(:)))+128;
imaga(msc);
offs=(n-n0)/(2*p.ds)+1;
hold on;

nClasses=max(d.rlnClassNumber);
for iCls=1:nClasses
    bMic=find(d.rlnClassNumber==iCls & p.iParts==micIndex);
    ourXs=d.rlnCoordinateX(bMic)/p.ds+offs(1);
    ourYs=d.rlnCoordinateY(bMic)/p.ds+offs(2);

%     bMic=bCls(iMics==micIndex);
    DrawBoxes(ourXs,ourYs,p.boxRad,p.colors(iCls,:));
end;
hold off;
[pa,nm,ex]=fileparts(name);
title([num2str(micIndex) '  ' nm ex '  ' num2str(p.micDefs(micIndex),3) 'um'],'interpreter','none');

end

function DrawBoxes(boxX,boxY,boxRadius,boxColor)

numB=numel(boxX);

sqx0=[-1 1  1 -1 -1 NaN]'*boxRadius;
sqy0=[-1 -1 1 1  -1 NaN]'*boxRadius;

nL=numel(sqx0);

boxesX=single(ones(nL,numB))*NaN;
boxesY=single(ones(nL,numB))*NaN;

for i=1:numB
        boxesX(:,i)=boxX(i)+sqx0;
        boxesY(:,i)=boxY(i)+sqy0;
end;
        boxesX=boxesX(:);
        boxesY=boxesY(:);
plot(boxesX(:),boxesY(:),'color',boxColor,'linewidth',1.5);

end
