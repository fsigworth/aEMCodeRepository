function mi=rspInsertAbsAmplitudes(mi,basePath,procPath)
% Given an mi file, load the corresponding rscc file and insert the
% absolute amplitude value from mxCCU

doDisplay=0;

if nargin<2
    basePath='';
end;
if nargin<3
    procPath=mi.procPath;
end;
nm=[AddSlash(basePath) AddSlash(procPath) mi.baseFilename 'rscc.mat'];
if ~exist(nm,'file')
    warning(['File doesn''t exist: ' nm]);
    return;
end;
load(nm);
if ~exist('mxCCU','var')
    error('Variable not loaded: mxCCU');
end;

ds=mi.imageSize(1)/size(mxCCU,1);
oldAmpsPresent=size(mi.particle.picks,2)>=9 && any(mi.particle.picks(:,9)>0);
nim=size(mi.particle.picks,1);
for i=1:nim
    coords=mi.particle.picks(i,:);
    if any(coords(3)==[16 32]) % optimized-manual or auto particle
        ix=round(coords(1)/ds+1);
        iy=round(coords(2)/ds+1);
        absAmp=mxCCU(ix,iy);
        if oldAmpsPresent
            oldStr=['  ** ' num2str(coords(9)) ' already set **'];
            disp([num2str(i,'%03d') '  ' num2str(coords(3)) '  ' num2str(absAmp) oldStr]);
        else
            mi.particle.picks(i,9)=absAmp;
            if doDisplay
                disp([num2str(i,'%03d') '  ' num2str(coords(3)) '  ' num2str(absAmp)]);
            end;
        end;
    end;
end;
