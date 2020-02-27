% rsReplaceMembraneModel
% Given an info structure mi, find and subtract vesicles
% and return the updated mi containing the vesicle coordinates.
% A vesicle model file is stored in Temp/ as <basename>v.mrc.
%
doRefineAnyway=1;  % even if previous refinement has been done.
useOkField=1;      % refine every vesicle for which ok is true.
doDownsampling=1;  % Downsample for speed
maxVesiclesToFit=inf;
disA=800;          % size of displayed/fitted window in angstroms
writeMiFile=1;     % Save the updated mi file
writeSubFile=1;    % Write a subtracted image into Temp/
setBasePath=1;     % Replace the mi.basePath with the path from the file selector.

% modelPath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Temp/';
% % modelName='sq02_10000ds2sym0all1bsc1binv0.mat';
% modelName='sq07_10003ds2sym0all1bsc0.5binv0alph0.015.mat';
modelPath='/Volumes/TetraData/EMWork/Hideki/121206p/AMPARliposome_antagonist_slot1/Temp/';
modelName='../sq07_10000ds4mbnModel.mat';
% modelName='sq07_10000ds4sym0all1bsc1binv0alph0.015.mat';
load([modelPath modelName]);

figure(1);
clf;
SetGrayscale;

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
mi.log=cell(0);
        mi0=mi;  % copy of original
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

        if setBasePath
            mi.basePath=rootPath;
            disp(['Reset mi.basePath to: ' mi.basePath]);
        end;
        
        mi0=mi;  % save the original, and load the new model
        disp(['loading the membrane model: ' modelName]);
        if ~exist([modelPath modelName],'file')
            [modelName modelPath]=uigetfile('*.mat','Select a membrane model');
            if numel(modelName)<1  %User clicked 'cancel'
                return
            end;
            modelPath=AddSlash(modelPath);
        end;
        load([modelPath modelName]);  % loads the vm structure
        if vm.pixA ~= mi.pixA
            disp(['Changing model pixel sizes: ' num2str([vm.pixA mi.pixA])]);
            mi.vesicleModel=meDownsampleVesicleModel(vm.vesicleModel,mi.pixA/vm.pixA);
        else
            mi.vesicleModel=vm.vesicleModel;  % insert the good model        
        end;
        
        %%
        
        % Get image and pixel sizes
        n=size(m,1);
        ds0=mi.imageSize(1)/n;  % downsampling factor of m
        pixA0=mi.pixA*ds0;    % pixel size of m
        ms=m;
        if doDownsampling
            % downsample to about 10A per pixel, yielding the image ms
            targetPixA=12;  % maximum pixel size
            ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
            if ns<n
                disp(['Downsampling to ' num2str(ns) ' pixels.']);
                ms=Downsample(m,ns);
            else
                ns=n;
            end;
            ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
            pixA=ds*mi.pixA;  %pixA in the image ms.
            ndis=NextNiceNumber(disA/pixA);
            
            %             nm=ceil(mi.imageSize(1)/ds;
        else  % use the original merged image scale
            ds=ds0;
            pixA=pixA0;
            ns=n;
            ndis=NextNiceNumber(disA/pixA)
        end;
        
        %%  Make model vesicles for the micrograph
        vs=meMakeModelVesicles(mi,ns);
% Least-squares fit the model + constant
        F=[vs(:) vs(:)*0+1];
        a=LinLeastSquares(F,ms(:));
        leastSquaresCoeffs=a
        vs1=a(1)*vs;
        mi0.vesicle.s=mi0.vesicle.s*a(1);  % change the amplitudes
        vs0=meMakeModelVesicles(mi0,ns);
        
        msub=ms-vs0;
        %%
        imacs(msub);
        drawnow;
         % Get the effective CTF from the merging.
%         H=ifftshift(meGetEffectiveCTF(mi,ns,ds));
        
        mi1=mi;  % this will be the output structure
        
% %%%%%%%% Actual fitting is done here
%         vfit=zeros(ns,ns);
%         figure(1)
%         disp('Fine fitting');
%         drawnow;
%         if ~isfield(mi.vesicle,'ok')
%             mi.vesicle.ok=1:numel(mi.vesicle.x);
%         end;
%         nVesicles=sum(mi.vesicle.ok>0);
%         nVesicles=min(nVesicles,maxVesiclesToFit);
%         nVesicles
%         for ind=1:nVesicles
%             if useOkField || mi1.vesicle(ind).ok
%                 figure(1);
%                 [mi1 diffIm vesFit]=rsFineFitVesicle(msub,mi1,mi0,ind,ndis,1,msub);
%                 subplot(2,2,2); title(fname{fileIndex},'interpreter','none');
% %                 figure(2); SetGrayscale;
% %                 [mi1 diffIm vesFit]=rsFineFitVesicleE(msub,mi1,mi0,ind,ndis,1,msub);
% %                 subplot(2,2,2); title(fname{fileIndex},'interpreter','none');
%                 
% % vindex=ind;
% % rsFitVesicleDeformation;
% % pause;
%             end;
%         end;
%         mi1.vesicle.refined=1;
%         
%%  Outputting
        if writeMiFile
            %%
            mi=mi1;
            outName=[infoPath mi.baseFilename 'mi.mat'];  % replace original
            if ~isfield(mi,'log')
                mi.log=cell(0);
            end;
            mi.log{numel(mi.log)+1}=['meRefineVesicleFits ' TimeStamp];
            save(outName,'mi');
            disp([outName ' saved']);
            
            %             Compute and store model vesicles
%             imacs(GaussFilt(m,.1));
            title('Original image');

            disp('Making the final vesicle models');
            vs1=meMakeModelVesicles(mi,n);

            imacs(GaussFilt(m-vs1,.1));
            title('Subtracted');
            drawnow;
            outVesName=[mi.tempPath mi.baseFilename 'v'];
            WriteMRC(vs1,pixA0,[outVesName '.mrc']);
            WriteJpeg(vs1,outVesName);
            imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
            disp([outVesName ' saved']);
        end;
        %               Store the subtracted micrograph
        if writeSubFile
            outSubName=[mi.tempPath mi.baseFilename 'mv'];
            WriteMRC(m-vs1,pixA0,[outSubName '.mrc']);
            WriteJpeg(m-vs1,outSubName);
%             imwrite(uint8(imscale(rot90(m-vs1),256,1e-3)),[outSubName '.jpg']);
            disp([outSubName ' saved']);
        end;
        
    else  % No vesicles have been found to refine
        if numel(mi.vesicle.x)<1
            disp('  ...no vesicles found');
        end;
        if isfield(mi.vesicle,'refined') && mi.vesicle.refined>0
            disp('  ...already refined.');
        end;
    end;
    disp(' ');
end; % for fileIndex
