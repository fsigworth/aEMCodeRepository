% meUpdateMiFiles12

cpe=16;
iCamera=1;  % F20 ultrascan
weights=[1 1 1];
writeFile=1;

% Ask for a set of original mi files.
[fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
if isnumeric(fname)
    return
end;
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);
nfiles=numel(fname);

%%
for findex=1:nfiles
    % Get the mi structure
    disp(['Reading ' infoPath fname{findex}]);
    load([infoPath fname{findex}]);
    if mi.version==11
        mi1=meCreateMicrographInfoStruct12;
        mi1.baseFilename=mi.baseFilename;
        if ~strcmp(mi.basePath,rootPath)
            disp(['changed basePath to: ' rootPath]);
        end;
        mi1.basePath    =AddSlash(rootPath);
        mi1.imagePath   =mi.imagePath;
        mi1.procPath    =mi.procPath;
        mi1.infoPath    =infoPath;  % use the present path
        if isfield(mi,'tempPath')
            mi1.tempPath    =mi.tempPath;
        end;
        mi1.imageFilenames=mi.imageFilenames;
        mi1.imageSize   =mi.imageSize;
        mi1.pixA        =mi.pixA;
        mi1.doses       =mi.doses(:)';  % row vector
        mi1.weights = mi1.doses*0+1;    % all ones
        mi1.kV          =mi.keV;
        mi1.camera      =  iCamera;
        mi1.cpe         =  cpe;
        mi1.ctf         =mi.ctf;
        mi1.mergeMatrix =mi.mergeMatrix;
        if isfield(mi,'mask')
            mi1.mask        =mi.mask;
        end;
        if numel(mi.noiseModelPars)>0 && ~any(isnan(mi.noiseModelPars))
            mi1.noiseModelCode   =mi.noiseModelCode;
            mi1.noiseModelPars =mi.noiseModelPars;
        end;
        mi1.vesicleModel=mi.vesicleModel;
        mi1.vesicle.x   =mi.vesicle.x(:);
        mi1.vesicle.y   =mi.vesicle.y(:);
        mi1.vesicle.r   =mi.vesicle.r(:);
        mi1.vesicle.s   =mi.vesicle.s(:);
        if isfield(mi.vesicle,'ok')
            nok=numel(mi.vesicle.ok);
            if nok>0
                mi1.vesicle.ok=repmat(mi.vesicle.ok(:),1,4);
            end;
        end;
        if isfield(mi.particle,'picks')
            mi1.particle.picks =mi.particle.picks;
        end;
        if isfield(mi.particle,'autopickPars')
            mi1.particle.autopickPars=mi.particle.autopickPars;
            mi1.particle.autopickPars(10)=0;  % extend to 10 elements
        end;
        if isfield(mi,'log');
            mi1.log=mi.log;
        end;
        if writeFile
            mi=mi1;
            disp(['Writing ' infoPath fname{findex}]);
            
            save([infoPath fname{findex}],'mi');
        end;
    else
        disp(['Version number is too high: ' num2str(mi.version)]);
    end;
    
end;

