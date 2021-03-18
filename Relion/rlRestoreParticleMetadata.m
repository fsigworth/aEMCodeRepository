% rlRestoreParticleMetadata

% Read a particles.star file A from the entire dataset (e.g. as used in
% Extraction) having the correct metadata. Read another star file B having selected particles.
% Make C, a copy of B, except+
%  - copy the optics group header from A to C
%  - copy metadata from A for each particle in B (recognized by same micrograph
% and coordinates) and place this into C. 

% The fields to be added or replaced are specified here:
% theFields={'rlnOpticsGroup'
%            'rlnGroupName'
%            'rlnGroupNumber'};
theFields={'vesicleRadius' ; 'vesiclePsi'};

% starAName='RSC/particleAll9_intens+frac_7505.star';
starAName='Class2D/job235/run_ct15_it020_data.star'
starBName='Class3D/job245/run_ct10_it025_data_original.star'
starCName='Class3D/job245/run_ct10_it025_data.star'

if ~exist(starAName,'file'), disp([starAName ' not found.']); end;
if ~exist(starBName,'file'), disp([starBName ' not found.']); end;
if ~exist(starCName,'file'), disp([starCName ' not found.']); end;
%%
disp(['Reading ' starAName]);
[nmA,aDat]=ReadStarFile(starAName);
aPts=aDat{2};
aNp=numel(aPts.rlnMicrographName);
%%
for i=1:numel(theFields)
    if ~isfield(aPts,theFields{i})
        disp(['This field is not found: ' theFields{i}]);
        return
    end;
end;
%%
disp(['Reading ' starBName]);
[nmB,bDat]=ReadStarFile(starBName);
bPts=bDat{2};
nbLines=numel(bPts.rlnMicrographName);
% Copy B into C
nmC=nmB;
cDat=bDat;
cPts=cDat{2};
    % Copy the optics info from A in t C
    cDat{1}=aDat{1};

    %% get the micrograph names from the subset file B
[bNames,~,bPartInds]=unique(bPts.rlnMicrographName);
bNumNames=numel(bNames);
disp([num2str(bNumNames) ' unique micrographs']);

% Get the full-set A micrograph names.
[aNames,~,aPartInds]=unique(aPts.rlnMicrographName);
disp([num2str(numel(aNames)) ' A micrograph refs']);

% Get iaLines, the A line numbers corresponding to each B line.
iaLines=zeros(bNumNames,1);
for i=1:bNumNames % look at each B micrograph name
    % find all the particles having that name
    ibParts=find(bPartInds==i);
    bXs=bPts.rlnCoordinateX(ibParts);
    bYs=bPts.rlnCoordinateY(ibParts);
    
    % find the correspoinding A name index
    j=find(strcmp(bNames{i},aNames),1);
    iaParts=find(aPartInds==j); % all the A particles with that name
    % now we find which coordinates match
    aXs=aPts.rlnCoordinateX(iaParts);
    aYs=aPts.rlnCoordinateY(iaParts);
    dists=hypot(bXs-aXs',bYs-aYs');
    [minDists,iaLocal]=min(dists,[],2);
    iaLines(ibParts)=iaParts(iaLocal);
 end;

 % At this point, iaLines(i) gives the line in A cooresponding to the ith line
 % in B.
 %% Now insert the values from A into B
nFields=numel(theFields);
for i=1:nFields
    fName=theFields{i};
    cPts.(fName)=aPts.(fName)(iaLines);
end;

cDat{2}=cPts;

% Write out the C star file.
disp(['Writing ' starCName]);
hdrText=['# version 30001; using metadata from ' starAName];
WriteStarFile(nmC,cDat,starCName,hdrText);
