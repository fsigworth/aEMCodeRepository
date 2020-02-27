% rsRefineVesicleFits
% Given an info structure mi, find and subtract vesicles
% and return the updated mi containing the vesicle coordinates.
% A vesicle model file is stored in Temp/ as <basename>v.mrc.
%
doRefineAnyway=1;  % even if previous refinement has been done.
useOkField=1;       % refine every vesicle for which ok is true.
showCoords=0;
doDownsampling=1;  % Downsample for speed
maxVesiclesToFit=inf;
disA=600;  % size of displayed window in angstroms
ndis=80;
writeFile=1;
writeSubFile=1;

modelPath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Temp/';
% modelName='sq02_10000ds2sym0all1bsc1binv0.mat';
modelName='sq02_10000ds4sym0all1bsc1binv0.mat';


% if nargin<1  % put up a file selector
[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
% pa
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);

for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    alreadyRefined=isfield(mi.vesicle,'refined') && mi.vesicle.refined>0
    %     check if there is something to do.
    if numel(mi.vesicle.x)>0 ...
        && (doRefineAnyway || (~alreadyRefined))

        iname=[mi.procPath mi.baseFilename 'm.mrc'];
        inamec=[mi.procPath mi.baseFilename 'mc.mrc'];
        if FileExists(iname)
            m=ReadEMFile(iname);
        elseif FileExists(inamec)
            m=ReadEMFile(inamec);
            iname=inamec;
        else
            [iname pa]=uigetfile('*m.mrc','Find the merged image');
            cd(pa);
            m=ReadEMFile(iname);
        end;
        
        %     Check that we have a temp directory
        if ~exist('Temp','dir')
            mkdir('Temp');
        end;
        mi.tempPath='Temp/';

%         % membrane model
%         if ~(isfield(mi,'vesicleModel') && numel(mi.vesicleModel)>1)
%             disp('Creating a generic membrane model');
%             % Don't have a valid model, make one.
%             vLipid=1.6;
%             thk=60;
%             rise=6;
%             % Create the model, which is sampled in units of the original pixel size.
%             nm0=ceil(30/mi.pixA)*2+1;  % array for vesicle model; 60A nominal
%             mi.vesicleModel=fuzzymask(nm0,1,thk/mi.pixA/2,rise/mi.pixA)*vLipid;
%             % units of V of inner potential
%         end;
% 
%         % check for old version of membrane scaling.
%         if max(mi.vesicleModel)>2
%             disp('Scaling down the vesicle model');
%             mi.vesicleModel=mi.vesicleModel/mi.pixA;
%         end;


% for i=1:numel(mi.ctf)
%     mi.ctf(i).B=mi.ctf(i).B/2;
% end;
% mi0.vesicle.s=mi0.vesicle.s*2;  % for some reason the GUI gives wrong values.

mi0=mi;  % save the original
disp(['loading ' modelName]);
load([modelPath modelName]);
mi.vesicleModel=vm.vesicleModel;  % insert the good model

        %%
        
        % Get image and pixel sizes
        n=size(m,1);
        ds0=mi.imageSize(1)/n;  % downsampling factor of m
        pixA0=mi.pixA*ds0;    % pixel size of m
        ms=m;
        if doDownsampling
            disp('Downsampling');
            % downsample to about 10A per pixel, yielding the image ms
            targetPixA=12;  % maximum pixel size
            ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
            if ns<n
                disp(['Downsampling to ' num2str(ns) ' pixels.']);
                ms=Downsample(m,ns);
            end;
            ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
            pixA=ds*mi.pixA;  %pixA in the image ms.
            ndis=NextNiceNumber(disA/pixA)

            %             nm=ceil(mi.imageSize(1)/ds;
        else  % use the original merged image scale
            ds=ds0;
            pixA=pixA0;
            ns=n;
            ndis=NextNiceNumber(disA/pixA)
        end;        
        
%         vm=Interpolate1(mi.vesicleModel,ns,1/ds)*pixA;  % downsampled vesicle model
%         mbnThickness=thk/pixA;
%         fms=fftn(ms);
        %%

        % get an initial subtraction
        vs=meMakeModelVesicles(mi0,ns);
        msub=ms-vs;
       %% 
        figure(1);
        SetGrayscale;
        subplot(2,3,1);
        imacs(msub);
        % Get the effective CTF from the merging.
        H=ifftshift(meGetEffectiveCTF(mi,ns,ds));
        
        
        % nv=numel(mi.vesicle.x);
        mi1=mi;  % this will be the output structure
%  rise=3;
%  ov=.2;
% thk=60;
% mi1=rsCreateVesicleModel(mi1,thk,rise,ov);

        indices=find(mi.vesicle.s>minS& mi.vesicle.s<maxS);
        nVesicles=numel(indices)
        mi1.vesicle.x=mi.vesicle.x(indices);
        mi1.vesicle.y=mi.vesicle.y(indices);
        mi1.vesicle.r=mi1.vesicle.r(indices);
        mi1.vesicle.s=mi1.vesicle.s(indices);
        mi1.vesicle.ok=ones(1,nVesicles);
        vfit=zeros(ns,ns);
        figure(1)
        disp('Fine fitting');
        drawnow;
        
nVesicles=min(nVesicles,maxVesiclesToFit);        
        
        for ind=1:nVesicles
            
            % un-subtract the vesicle
            %     v=meMakeModelVesicles(mi,ns,ind);
            %     msub2=msub+v;
            %     toc
            % Do the fine fitting
%              mi1=rsFineFitVesicle2(msub,mi1,ind,ndis,1,msub);  % take default ndis=128.
            mi1=rsFineFitVesicle(msub,mi1,mi0,ind,ndis,1,msub);  % take default ndis=128.
            mi1.vesicle.ok(ind)=1;
%             subplot(2,2,1); title(ind);
            subplot(2,2,2); title(fname,'interpreter','none');
            %     vfit=vfit+meMakeModelVesicles(mi1,ns,ind);
            %     toc
            %     figure(2);
            %     subplot(1,1,1);
            %     SetGrayscale;
            %     imacs(ms-vfit);
            %     title(ind);
            %     drawnow;
            if showCoords
                disp([ind mi1.vesicle.x(ind) mi1.vesicle.y(ind) mi1.vesicle.r(ind) mi1.vesicle.s(ind)]);
            end;
        end;
        mi1.vesicle.refined=1;
        
        %%
        if writeFile
            %%
            mi=mi1;
            %     outName=[infoPath mi.baseFilename 'smi.mat'];
            outName=[infoPath mi.baseFilename 'mi.mat'];  % replace original
            if ~isfield(mi,'log')
                mi.log=cell(0);
            end;
            mi.log{numel(mi.log)+1}='meRefineVesicleFits';
            save(outName,'mi');
            disp([outName ' saved']);
            vs1=meMakeModelVesicles(mi,n);
  %%
  outVesName=[mi.tempPath mi.baseFilename 'v41'];
            WriteMRC(vs1,pixA0,[outVesName '.mrc']);
            imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
            disp([outVesName ' saved']);
            if writeSubFile
                outSubName=[mi.tempPath mi.baseFilename 'mv41'];
                WriteMRC(m-vs1,pixA0,[outSubName '.mrc']);
                imwrite(uint8(imscale(rot90(m-vs1),256,1e-3)),[outSubName '.jpg']);
                disp([outSubName ' saved']);
            end;
            
        end;
    else  % Nothing to do.
        if numel(mi.vesicle.x)<1
            disp('  ...no vesicles found');
        end;
        if isfield(mi.vesicle,'refined') && mi.vesicle.refined>0
            disp('  ...already refined.');
        end;
    end;
    
end; % fileIndex