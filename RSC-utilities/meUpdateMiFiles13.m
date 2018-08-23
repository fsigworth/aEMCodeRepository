% meUpdateMiFiles13
% update the mi structure to the version 13 standards.
% This can involve re-computing the merged images to fix the scaling.

cpe=16
ds=2;  % merged iaage downsampling
defaultKV=200;
iCamera=1;  % F20 ultrascan
weights=[1 1 1]
writeMiFile=1;
forceRemerge=1  % force all images to be re-merged.

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
disp('Updating mi files to version 13.');
for findex=1:nfiles
    % Get the mi structure
    disp(['Reading ' infoPath fname{findex}]);
    load([infoPath fname{findex}]);
    disp(['  mi version ' num2str(mi.version)]);
    mi.basePath=rootPath;
    %     Get the merged image
    disp(['Reading ' mi.procPath mi.baseFilename 'm.mrc']);
        m=ReadMRC([mi.procPath mi.baseFilename 'm.mrc']);
        autoRemerge=false;
mi.imagePath='Micrograph/';

        switch mi.version
            case {11,12}
                mi1=meCreateMicrographInfoStruct13;
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
                if isfield(mi,'keV') && mi.keV>0
                    mi1.kV=mi.keV;
                elseif isfield(mi,'kV') && mi.kV>0
                    mi1.kV=mi.kV;
                else
                    mi1.kV=defaultKV;
                end;
                
                mi1.camera      =  iCamera;
                if isfield(mi,'cpe')
                    mi1.cpe          =mi.cpe;
                else
                    mi1.cpe = cpe;
                end;
                
                mi1.ctf         =mi.ctf;
                if isfield(mi1.ctf,'pixA')
                    mi1.ctf=rmfield(mi1.ctf,'pixA');
                end;
                if isfield(mi1.ctf,'res')
                    mi1.ctf=rmfield(mi1.ctf,'res');
                end;
                
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
                    noks=size(mi.vesicle.ok);
                    if noks(1)>0 && noks(2)<4
                        mi1.vesicle.ok=repmat(mi.vesicle.ok(:,1),1,4);
                    elseif noks(1)<1
                        mi1.vesicle.ok=true(numel(mi1.vesicle.x),4);
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
                
                mi=mi1;
                if writeMiFile
                    disp(['Writing ' infoPath fname{findex}]);
                    save([infoPath fname{findex}],'mi');
                end;
                if numel(mi.vesicle.s>0) && sum(~isnan(mi.vesicle.s))>0
                    medianVesAmplitude=median(mi.vesicle.s)
                    autoRemerge= (medianVesAmplitude > .008); % probably old scaling
                end;
            case 13
                disp('mi file is already version 13.');
                if numel(mi.vesicle.s>0) && sum(~isnan(mi.vesicle.s)>0)
                    medianVesAmplitude=median(mi.vesicle.s(~isnan(mi.vesicle.s)))
                    autoRemerge= (medianVesAmplitude > .008); % probably old scaling
                end;
                if writeMiFile
                    disp(['Writing ' infoPath fname{findex}]);
                    save([infoPath fname{findex}],'mi');
                end;

            otherwise
                disp(['Version number is out of range: ' num2str(mi.version)]);
        end;
        if autoRemerge || forceRemerge
            disp('Re-merging the image');
            m=sum(meMakeMergeImageSet(mi,mi.cpe,ds),3);  % create a new merged image
            procname=[mi.procPath mi.baseFilename];  % Use the processed image directory
            WriteMRC(m,ds*mi.pixA,[procname 'm.mrc']);
            figure(1); clf; SetGrayscale;
            imacs(m);
            title(procname,'interpreter','none');
            drawnow;
            disp(['Wrote the merged image: ' procname 'm.mrc']);
            nd=round(mi.imageSize/8);
            vd=meMakeModelVesicles(mi,nd);
            md=Downsample(m,nd);
            vesScale=(md(:)'*vd(:))/(vd(:)'*vd(:));
            disp(['Vesicle scale: ' num2str(vesScale)]);
            if vesScale>1.5 || vesScale <0.7
                mi.vesicle.s=mi.vesicle.s*vesScale;
                disp(['Updating ' infoPath fname{findex}]);
                save([infoPath fname{findex}],'mi');                
            end;
        end;
        
        disp(' ');
end;

