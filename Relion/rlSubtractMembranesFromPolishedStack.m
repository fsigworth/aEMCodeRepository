% rlSubtractMembranesFromPolishedStack
% We read the star file from the last Refine3D, and from the subsequent
% polishing step. We'll call these Ref and Pol.
% - In Ref.star we get the .mrcs data files for the rlnImageName and
% rlnReconstructImageName.
% - We subtract pairs of images in thes files, to get the membrane densities.
% - We get the rlnImageName from Pol.star, and read it from the *shiny.mrcs file. We
%  subtract the membrane density and write out the result into a new *shiny_v.mrcs
%  into the same Polish/job directory.
% - We create a new shiny_v+unsub.star file which points to the
% new *shiny_v.mrcs for rlnImageName, and to the old *shiny.mrcs for the rln
% ReconstructImageName. We'll write this back into the Polish/job directory. 
%
% In this simple version we assume a 1-to-1 correspondence between the
% lines in the Ref and Pol files.

refStarName='Refine3D/job312/run_data.star'; % sub + unsub stack
polishDir='Polish/job318/'; % where to find the input shiny.star
doWrite=1; % Write the output files back into polishDir!

displayCount=10; % 1 in this many files get displayed.
% The following should always be right for the Polish/ directory
polStarName=[polishDir 'shiny.star'];
outStarName=[polishDir 'shiny_v+unsub.star'];

% Get the star files. To save time, we first check for the presence of the Pol file.
if ~exist(polStarName,'file')
    disp(['The input file ' polStarName ' can''t be found.']);
    return
end;

disp(['Reading Refine3D reference: ' refStarName]);
[refNm,refDat]=ReadStarFile(refStarName);
refD=refDat{2}; % 2nd block

disp(['Reading the polished star file: ' polStarName]);
[polNm,polDat]=ReadStarFile(polStarName);
vPolDat=polDat; % copy for output file
vPolD=vPolDat{2}; % copy the data struct
nPolRows=numel(vPolD.rlnImageName);


disp('Decoding ref names');
[vInds,vFileNames]=rlDecodeImageName(refD.rlnImageName);
[uInds,uFileNames]=rlDecodeImageName(refD.rlnReconstructImageName);

%% --get the unique stack file names---
[vRefNames,vFirstRows,vNameInds]=unique(vFileNames,'stable'); % leave them in original order
uRefNames=uFileNames(vFirstRows);
numUnique=numel(vFirstRows);

% The Pol imageNames are the unsubtracted ones.
% Immediately insert the rlnReconstructImageName, which points to the
% present unsub particle image.
vPolD.rlnReconstructImageName=vPolD.rlnImageName;
[~,pUFileNames]=rlDecodeImageName(vPolD.rlnImageName(vFirstRows));

nVRows=numel(refD.rlnImageName);

if (nPolRows ~= nVRows)
    disp(['Inconsistent num Ref and Polish rows: ' num2str([nVRows nPolRows])]);
    return
end;

% Construct the subtracted Pol names and make the files.
pVFileNames=cell(numUnique,1);
disp([num2str(numUnique) ' files to write.']);
for i=1:numUnique
    inUName=pUFileNames{i};
    [pa,nm,ex]=fileparts(inUName);
    outVName=[AddSlash(pa) nm '_v' ex]; % add _v to the filename.
    pVFileNames{i}=outVName;
    [uPolParts,s]=ReadMRC(inUName);
    uRefParts=ReadMRC(uRefNames{i});
    vRefParts=ReadMRC(vRefNames{i});
    mModels=uRefParts-vRefParts; % subtract all the particles
    vPolParts=uPolParts-mModels;
    if mod(i,displayCount)==1 % we'll display one particle's results.
        n=size(uPolParts,1);
        disDs=NextNiceNumber(ceil(n/64),2); % a power of 2
        nDis=NextNiceNumber(n/disDs);
        subplot(221);
        imags(Downsample(uPolParts(:,:,1),nDis));
        title([nm '_v' ex],'interpreter','none');
        subplot(222);
        imags(Downsample(mModels(:,:,1),nDis));
        subplot(223);
        imags(Downsample(vPolParts(:,:,1),nDis));
        title('Subtracted');
        drawnow;
    end;
    if doWrite
        WriteMRC(vPolParts,s.pixA,outVName);
        disp([num2str(i) '  ' outVName]);
    end;
end;

for i=1:nPolRows
    stackName=pVFileNames{vNameInds(i)};
    vPolD.rlnImageName{i}=[sprintf('%06d',vInds(i)) '@' stackName];
    vPolD.rlnMicrographName{i}=refD.rlnMicrographName{i};
end;

vPolDat{2}=vPolD;
if doWrite
    WriteStarFile(polNm,vPolDat,outStarName);
    disp([outStarName ' written.']);
end;


