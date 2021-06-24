% rlShowParticlePicks.m
% Based on a particles star file, load a micrograph and show the particle
% locations.
% Based roughly on DemoPicker

starPath='RSC9/';
starName='particles_10k_u.star';

% micPath='Merged_sms/';
% micSuffix='ms.mrc';
% micTrim=2; % no. of chars to trim from micrograph name (not extension)
% ds=4;
% disDs=1;

micPath='MotionCorr/job007/Movies/';
micSuffix='.mrc';
micTrim=2; % no. of chars to trim from micrograph name (not extension)
ds=1;    % Existing downsampling of the micrographs
disDs=4; % Downsampling here for display

disSD=1e-3; % Display clips this fraction of the intensity range.
disFc=.1; % Display Gaussian filter, in A^-1
lfAmp=.4; % CTF compensation fraction

lineWidth=2;

[nms,das,ok]=ReadStarFile([starPath starName]);
d=das{end};
np=numel(d.rlnMicrographName);
[uNames,firstInds,nameInds]=unique(d.rlnMicrographName,'stable');
nMics=numel(uNames);
disp([num2str(nMics) ' micrographs.']);

ind=1;
doLoad=1;
boxSizeA=150;
figure(1);
ax=gca;

% For displaying multiple classes, the colors are indexed by class group
colors= [0 1 .3;
         .5 1 0];
b='.';
while b~='q'
    doLoad=1; % by default, load a new micrograph
    switch b
        case {'n' ' '} % Get next micrograph
            ind=ind+1;
            if ind>nMics
                beep;
                ind=nMics;
            end;
        case {'N' 'p'} % Get the previous micrograph
            ind=ind-1;
            if ind<1
                beep;
                ind=1;
            end;
        case 'g' % go to file index
            ind=MyInput('go to file index ',ind);
            ind=min(nMics,max(1,ind));
%         case 'b' % change box display
%             boxesOn=boxesOn+1;
%             if boxesOn>2
%                 boxesOn=0;
%             end;
%             doLoad=0;
    end;
    if doLoad
        iLine=firstInds(ind);
        [ok,~,mi]=rlStarLineToMi(nms,das,iLine);
        [pa,nm]=fileparts(uNames {ind});
        micFilename=[micPath nm(1:end-micTrim) micSuffix];
        str=[num2str(ind) '  ' micFilename];
        disp(str);
        [m0,s]=ReadMRC(micFilename);
        n0=size(m0);
        n1=NextNiceNumber(n0,5,8);
        if any(n1~=n0) % need to pad the image with its mean.
            me=mean(m0(:));
            m0=m0-me;
            m0(n1(1),n1(2))=0; % expand the image.
            m0=m0+me;
        end;
        mi.padImageSize=ds*n1;
        if ind==1
            padImageSize=mi.padImageSize
        end;
        pixA=s.pixA*disDs;
        m1=GaussFilt(Downsample(m0,n1/disDs),disFc*pixA);
        if lfAmp>0
            m1=meCTFInverseFilter(m1,mi,lfAmp);
        end;
    end;
    
    % Make the display
    mImg=imscale(m1,256,disSD);
    imaga(mImg);
    title(str,'interpreter','none');
    
    pInds=find(nameInds==ind);
    nParts=numel(pInds);
    if nParts>0 % We have boxes to draw
        disX=d.rlnCoordinateX(pInds)/(ds*disDs)+1;
        disY=d.rlnCoordinateY(pInds)/(ds*disDs)+1;
        class=ones(nParts,1,'single');
        boxSize=boxSizeA/(ds*disDs);
        boxXs=[-1 1 1 -1 -1]*boxSize/2;
        boxYs=[-1 -1 1 1 -1]*boxSize/2;
        hold on;
        for i=1:nParts
            plot(boxXs+disX(i),boxYs+disY(i),'-','color',colors(class(i),:),'linewidth',lineWidth);
        end;
        hold off;
        drawnow;
        [~,~,b]=Myginput(1);
    end;
end;
    disp('done.');
