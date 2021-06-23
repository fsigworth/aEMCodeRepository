% rlShowParticlePicks.m

starPath='';
starName='';
micPath='';
micSuffix='ms.mrc';
disSD=1e-3;

ds=4;
disDs=1;
disFc=.1; % in A^-1

[nms,das,ok]=ReadStarFile([starPath starName]);
d=das{end};

[uNames,~,nameInds]=unique(d.rlnMicrographName,'stable');
nMics=numel(uNames);
disp([num2str(nMics) ' micrographs.']);

ind=1;
doLoad=1;
boxSizeA=200;
figure(1);
ax=gca;
colors=ax.ColorOrder; % handles 7 colors.

while b~='q'
    switch b
        case {'n' ' ' 'R'} % Get next file or repeat
            ind=ind+1;
            if ind>nl
                beep;
                ind=nl;
            end;
            doLoad=1;
        case {'N' 'p'} % Get the previous file
            ind=ind-1;
            if ind<1
                beep;
                ind=1;
            end;
            doLoad=1;
        case 'g' % go to file index
            ind=MyInput('go to file index ',ind);
            ind=max(1,ind);
            doLoad=1;
        case 'b' % change box display
            boxesOn=boxesOn+1;
            if boxesOn>2
                boxesOn=0;
            end;
 
            
if doLoad
    [pa,nm]=fileparts(rlnMicrographName{ind});
    micFilename=[micPath nm micSuffix];
    [m0,s]=ReadMRC(micFilename);
    n0=size(m0);
    n1=NextNiceNumber(n0,5,8);
    if any(n1~=n0) % need to pad the image with its mean.
        me=mean(m0(:));
        m0=m0-me;
        m0(n1(1),n1(2))=0; % expand the image.
        m0=m0+me;
    end;
    pixA=s.pixA*disDs;
    m1=GaussFilt(Downsample(m0,n1/disDs),disFc*pixA);
end;

% Make the display
mImg=imscale(m1,256,disSD);
imaga(mImg);

pInds=find(nameInds==ind);
nParts=numel(pInds);
if nParts>0 % We have boxes to draw
disX=d.rlnCoordinateX(pInds)/(ds*disDs)+1;
disY=d.rlnCoordinateY(pInds)/(ds*disDs)+1;
class=ones(nParts,1,'single');
boxSize=boxSizeA/(ds*disDs);
boxXs=[-1 1 1 -1]*boxSize/2;
boxYs=[-1 -1 1 1]*boxSize/2;
hold on;
for i=1:nParts
    plot(boxXs+disX(i),boxYs+disY(i),'-','color',color(class(i),:));
end;
hold off;

[~,~,b]=Myginput(1);
end;

