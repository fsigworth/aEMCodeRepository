% rlAddReconstructImageName.m
% add the column to a particles.star file.
% We read the particles.star files from the two Extract jobs created from
% our newly-generated particles_u and particles_v files (from 
% rlMisToParticleStar). We depend on there being a line-for-line
% correspondence between the rlnImageNames in the two files.
% We match filenames in the input file (which comes e.g. from selection or
% refinement) with those in the _v star file. Then we take the
% corresponding filename from the _u file and construct the complete unsub  
% image names and place them into the new rlnReconstructImageName field of
% the input file. The modified input file is written out with "+unsub"
% appended to the name. Given an input file run_data.star, the output will
% be run_data+unsub.star

% particles.star name
% inStarName='Select/job023/particles.star';
inStarName='RSC1_C25-3/particles_iso_all.star';
refVStarName='Extract/job012/particles.star';
refUStarName='Extract/job013/particles.star';

insertReconstructImage=1; % Yes, do put in the field
doSwap=0; % Put the _u particles into rlnImageName, etc?
[pa, nm, ex]=fileparts(inStarName);
outStarName=[pa filesep nm '+unsub' ex];

% % out of desparation I"m trying swapping the imageName and reconstrImageName
% % fields:
% outStarName=[pa filesep nm '+uswapd' ex]; doSwap=1;

% Get the pair of micrograph.star files. We assume that the unsub and
% subtracted micrograph names are on the same rows in the files.
if ~exist(inStarName,'file')
    disp(['The input file ' inStarName ' can''t be found.']);
    return
end;

disp(['Reading V reference: ' refVStarName]);
[nmv,datv]=ReadStarFile(refVStarName);
vImgNames=datv{2}.rlnImageName;

disp(['Reading U reference: ' refUStarName]);
[nmu,datu]=ReadStarFile(refUStarName);
uImgNames=datu{2}.rlnImageName;

nUNames=numel(uImgNames);
if numel(vImgNames)~=nUNames
    disp('Numbers of rows don''t match. Exiting.');
    return
end;

disp('Decoding refVStar names.')
[~,vFileNames]=rlDecodeImageName(vImgNames);
[vFileUniques,vFirst,vLookup]=unique(vFileNames);

disp('Decoding refUStar names.')
[~,uFileNames]=rlDecodeImageName(uImgNames);

% Get the input particle.star file.
disp(['Reading the input file: ' inStarName]);
[nmp,datp]=ReadStarFile(inStarName);
nparts=numel(datp{2}.rlnImageName);
dOut=datp{2}; % Copy the particle data.

disp(['Decoding ' num2str(nparts) ' particle names']);
[pSlices,inFileNames]=rlDecodeImageName(dOut.rlnImageName);
[pFileUniques,~,pLookup]=unique(inFileNames);

% match the particle filenames with refVStar
nu=numel(pFileUniques);

matchFileInds=zeros(nu,1);
for i=1:numel(pFileUniques)
    q=strcmp(pFileUniques{i},vFileUniques);
    ind=find(q);
    if numel(ind)>1
        disp(['Duplicated micrograph name? ' pFileUniques{i}]);
    elseif numel(ind)<1
        disp(['No matching micrograph name: ' pFileUniques{i}]);
        continue;
    end;
    matchFileInds(i)=ind;
end;

%% for each pLookup, get matchFileInds(pLookup) which point into
% vFileUniques. vFirst(matchFileInds(pLookup)) points to the uFileName.

dOut.rlnReconstructImageName=cell(nparts,1);
disp('Creating the new particle names');
for i=1:nparts
    newNameInd=vFirst(matchFileInds(pLookup(i)));
    newName=[num2str(pSlices(i),'%06d') '@' uFileNames{newNameInd}];
    dOut.rlnReconstructImageName{i}=newName;
end;

if doSwap
%%swap the fields
temp=dOut.rlnReconstructImageName
% Include this field only if we want it.
if insertReconstructImage
    dOut.rlnReconstructImageName=dOut.rlnImageName;
else
    dOut=rmfield(dOut,'rlnReconstructImageName');
end;
dOut.rlnImageName=temp;
%%%%
end;

disp(['Writing the output file ' outStarName]);
dato=datp;
dato{2}=dOut;
WriteStarFile(nmp,dato,outStarName);
disp('done.')



% Primitive old code.

% 
% 
% 
% inStarName='Select/job191/particles.star';
% outStarName='Select/job191/particles+unsub.star';
% 
% disp(['Reading ' inStarName]);
% imgNames=dat{2}.rlnImageName;
% 
% whos imgNames
% %%
% outNames=imgNames;
% nim=numel(imgNames);
% for i=1:nim
%     [partInd,micName]=rlDecodeImageName(imgNames{i});
%     [pa,nm,ex]=fileparts(micName);
%     nm(end)='u'; % swap the final character from 'v' to 'u'
%     name=imgNames{i};
%     if strndcmp(name(end-4:end),'.mrcs')
%         name(end-5)='u'; % mark unsubtracted
%         disp(name);
%         return
%     end;
% end;


% Simple old code...
%
% % rlAddReconstructImageName
% % add the column to a particles.star file
% 
% inStarName='Select/job191/particles.star';
% outStarName='Select/job191/particles+unsub.star';
% 
% disp(['Reading ' inStarName]);
% [nm,dat]=ReadStarFile(inStarName);
% imgNames=dat{2}.rlnImageName;
% 
% whos imgNames
% %%
% outNames=imgNames;
% nim=numel(imgNames);
% for i=1:nim
%     [partInd,micName]=rlDecodeImageName(imgNames{i});
%     [pa,nm,ex]=fileparts(micName);
%     nm(end)='u'; % swap the final character from 'v' to 'u'
%     name=imgNames{i};
%     if strndcmp(name(end-4:end),'.mrcs')
%         name(end-5)='u'; % mark unsubtracted
%         disp(name);
%         return
%     end;
% end;