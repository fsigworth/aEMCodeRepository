% rtOnlineDisplayNoGraphics.m
inWorkingDir=0;  % not finished
tempDir='temp/';
imageDir='Micrograph/';
jpegDir='Jpeg/';
infoDir='Info/';
procDir='Merged/';
ourRefName='CountRefLocal.mrc';
gainRefRot=0;
kV=200;
searchMode=2; % 1: operate on latest; 2 start from beginning. 
pixA=1.247;
cpe=20;
times=0;
defoci=0;
maxRepeat=100; % 500 seconds max. waiting before quitting.
maxRepeat=1;

mc2Pars=struct;
mc2Pars.throw=1;   % for MC2: default is zero.
mc2Pars.trunc=0;  % default is zero
mc2Pars.gainRefRot=0;

disp('Getting the first movie');
[firstName,inPath]=uigetfile('*.*','Select a movie file');
if ~inWorkingDir
    cd(inPath);  % Might need to use a base path
    inPath='';
end;
[pa,nm,ex]=fileparts(firstName);
namePattern=ex;

%%
% Look for a gain reference
d=dir(inPath); % same as movie dir
found=0;
for i=3:numel(d)
    [~,rnm,rex]=fileparts(d(i).name);
    if strcmp(rex,'.dm4') && numel(strfind(lower(rnm),'ref'))
        refName=d(i).name;
        refPath=GetRelativePath(inPath);
        found=i;
        break;
    end;
end;
if found
    gainRefName=[refPath refName];
else
    disp('Getting the refernce');
    [refName,refPath]=uigetfile('*.dm4','Select the reference');
    if isnumeric(refName) % nothing selected
        gainRefName='';
    else
        gainRefName=[GetRelativePath(refPath) refName];
    end;
end;
if numel(gainRefName)>1
    ref=ReadEMFile(gainRefName);
    WriteMRC(ref,0,ourRefName);
end;

if numel(gainRefName)>1
disp(['Gain reference: ' gainRefName ' copied to ' ourRefName]);
else
    disp('No gain reference');
    ourRefName='';
end;

nameList=cell(0);
timeList=[];

CheckAndMakeDir(tempDir, 1);
CheckAndMakeDir(imageDir, 1);
CheckAndMakeDir(jpegDir, 1);
CheckAndMakeDir(infoDir,1);
CheckAndMakeDir(procDir,1);

mi0=meCreateMicrographInfoStruct14;
mi0.pixA=pixA;
% mi0.movieFilename{1}=nameList{end};
mi0.tempPath=tempDir;
mi0.imagePath=imageDir;
mi0.procPath=procDir;
mi0.jpegPath=jpegDir;
mi0.gainRefName=ourRefName;
mi0.gainRefRot=gainRefRot;
mi0.kV=kV;
mi0.cpe=cpe;
mi0.camera=5;
mi0.moviePath='';
    if mi0.kV>200
        mi0.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi0.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;

%% ---------------------Main loop---------------
doRepeat=1;
miIndex=0;

while doRepeat
    ok=0;
    repeatCount=0;
    while ~ok && repeatCount<maxRepeat
        [nameList,timeList,ok]=rtFindNextMovie(nameList,timeList,namePattern,searchMode);
        if ~ok
            sprintf('.');
            pause(5);
        end;
        repeatCount=repeatCount+1;
    end;
    doRepeat=repeatCount<=maxRepeat;
    
    miIndex=miIndex+1;
    times(miIndex)=timeList(end);
    

    mi=mi0;
    mi.movieFilename{1}=nameList{end};
    [pa,nm,ex]=fileparts(mi.movieFilename{1});
    mi.baseFilename=nm;
    mi.imageFilenames{1}=[nm 'ali.mrc'];
    
    % ------Run MC2 here------
    mi=rtMC2Runner(mi,mc2Pars);
    if isnumeric(mi) % MC2 run was unsuccessful.
        disp('MotionCor2 run was unsuccessful.');
%        return;
        break;
    end;
    m=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    n0=size(m);
    n=NextNiceNumber(n0,5,8);
    m1=Crop(m,n,0,mean(m(:))); % expand to nice size
    m1s=Downsample(m1,n/4);
    
%%%%%
    s=struct;
    s.exec='mysubplot(122)';
    s.image=m1s;
    s.axis='off';
    s.title=mi.movieFilename{1};
    disDat.mainImage=s;
    
%     mysubplot(122);
%     imags(m1s);
%     axis off;
%     title(mi.movieFilename{1}, 'interpreter','none');
%%%%%
    s=struct;
    s.exec='mysubplot(2,4,1)';
    s.ploty=mi.frameShifts{1};
    s.exec('grid on');
    s.ylabel('Shift, pixels');
    disDat.frameShifts=s;
%     mysubplot(2,4,1);
%     plot(mi.frameShifts{1});
%     grid on;
%     ylabel('Shift, pixels');
%     drawnow;
    
    gPars=struct;
    [mi,epaVals,ctfImage,ctfVals]=rtGctfRunner(mi,gPars);
    if isnumeric(mi)
        disp('Skipping');
        break
    else
        disp(['defocus ' num2str(mi.ctf(1).defocus)]);
    end;
    
    % put this onto the drift plot  
    %%%%
    disDat.frameShifts.title=['Res limit: ' num2str(ctfVals.RES_LIMIT,3) 'Å'];
%     title(['Res limit: ' num2str(ctfVals.RES_LIMIT,3) 'Å']);

    %%%%%
    s=struct;
    s.exec='mysubplot(242)';
    s.image=ctfImage;
    s.axis='off';
    s.title=['\delta = ' num2str(mi.ctf.defocus,3) '\mum ' ...
        ' astig = ' num2str(abs(mi.ctf.deltadef),3) '\mum'];
    disDat.ctfImage=s;
%     mysubplot(2,4,2);
%     imags(ctfImage);
%     axis off;
%     title(['\delta = ' num2str(mi.ctf.defocus,3) '\mum  astig = ' num2str(abs(mi.ctf.deltadef),3) '\mum']);
%     
    %%%%%
    s=struct;
    s.exec='mysubplot(425)';
    mysubplot(4,2,5);
    scl4=1/max(abs(epaVals.epaBkgSub));
    s.xplot=1./epaVals.resolution;
    s.yplot=[epaVals.ctfSim.^2 scl4*epaVals.epaBkgSub epaVals.ccc 0*epaVals.ccc];
    disDat.epa=s;
    
%     plot(1./epaVals.resolution,[epaVals.ctfSim.^2 scl4*epaVals.epaBkgSub epaVals.ccc 0*epaVals.ccc]);
    
    defoci(miIndex)=mi.ctf.defocus;
    %%%%%
    s=struct;
    s.exec='subplot(427)';
%     mysubplot(4,2,7);
    relTimes=(times-min(times))*24;
    s.xplot=relTimes;
    s.yplot=defoci;
%     plot(relTimes,defoci);
    s.ylabel='Defocus \mum';
    s.xlabel='Time, hr';
    
%     ylabel('Defocus, \mum');
%     xlabel('Time, hr');
    disDat.defocusPlot=s;
    
%     drawnow;
    
    %%%%
    jpegName=[mi.jpegPath mi.baseFilename '.jpg'];
    disDat.print.eval=['print(' jpegName '''-djpeg'')'];
%     print(jpegName,'-djpeg');
    save([mi.jpegPath mi.baseFilename 'dis.dat']);

    WriteMiFile(mi,[infoDir mi.baseFilename 'mi.txt']);
    
    
end;  % big while



function [nameList,timeList,ok]=rtFindNextMovie(nameList,timeList,pattern,mode)
ok=0;
d=dir;
j=0;
names=cell(0);
times=[];
for i=1:numel(d)
    nameOk=numel(strfind(d(i).name,pattern));
    if ~d(i).isdir && nameOk
        j=j+1;
        names{j}=d(i).name;
        times(j)=d(i).datenum;
    end;
end;
if j==0
    return
end;
switch mode
    case 1  % just pick the last
        [newTime,newInd]=max(times);
        if numel(timeList)<1 || newTime>max(timeList)
            timeList(end+1,1)=newTime;
            nameList{end+1,1}=d(newInd).name;
            ok=1;
        end;
    case 2  % starting from the beginning, work forward.
        if numel(timeList)>0
            oldTime=timeList(end);
        else
            oldTime=0;
        end;
        [sortedTimes,sortedIndices]=sort(times);
        newInd=sortedIndices(find(sortedTimes>oldTime,1,'first'));
        if numel(newInd)>0
            timeList(end+1,1)=times(newInd);
            nameList{end+1,1}=names{newInd};
            ok=1;
        end;
end
end