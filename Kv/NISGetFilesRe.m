function [names,pd,ptrsI,ligandLabels,sites,P]=NISGetFilesRe;

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia/')
% names.map='211021/cryosparc_P1_J340_005_volume_map_sharp.mrc';
% names.pdb='211021/I_10-18-coot-21_real_space_refined_003 (1).pdb';
% names.outDir='220216/';
% names.dataName='NISMapDataI.mat'; % output file.

names.map='211116/cryosparc_P8_J137_003_REO_map_sharp.mrc'; % The ReO map
names.pdb='211116/REO_11-16-21-coot-9_real_space_refined_009.pdb';
names.dataName='NISMapDataREO.mat'; % output file.
names.outDir='220216/';

P=[67 2.3];  % rotation, z-translation
% Get the displacement of Na2 from Na1 from the I pdb.
NaDiff=[169.415 166.937 154.764] - [169.814 168.578 158.465]

disp(['Loading the pdb file 5: ' names.pdb]);
pd=pdb2mat(names.pdb);

%% locate the ions
disp('Locating the ions in pdb 5');

NaPtrs=find(strcmp('Na',pd.element));
if numel(NaPtrs)<2 % We'll put in the location of Na2 from the I model
    % In line 2. The 
    NaPtrs(2)=NaPtrs(1)+1; % fake pointer for Na2
    
    pd.element{NaPtrs(2)}='Na';
    pd.atomNum(NaPtrs(2))=NaPtrs(2);
    pd.chainID{NaPtrs(2)}=pd.chainID{NaPtrs(1)};
    pd.resNum(NaPtrs(2))=pd.resNum(NaPtrs(1))+1;

    pd.X(NaPtrs(2))=pd.X(NaPtrs(1))+NaDiff(1);
    pd.Y(NaPtrs(2))=pd.Y(NaPtrs(1))+NaDiff(2);
    pd.Z(NaPtrs(2))=pd.Y(NaPtrs(1))+NaDiff(3);
end;
IPtrs=find(strcmp('Re',pd.element));
QPtrs=find((pd.resNum==72 | pd.resNum==94) & strcmp('CA',pd.atomName));
disp(['  Gln pointers ' num2str(QPtrs)]);

ptrsI=[IPtrs NaPtrs QPtrs];
nL=numel(ptrsI);

ligandLabels=cell(nL,1);
for i=1:nL
    txt=[pd.element{ptrsI(i)} ' ' num2str(pd.atomNum(ptrsI(i)))];
    disp(txt);
    ligandLabels{i}=txt;
end;

% Find the binding site atoms


%%
% Look up the coordinating atoms in the binding sites.
% We return the pdb file atom nos. (atomLines) of the closest atom in the
% residue.

nIons=1+numel(NaPtrs);
sites=struct;
for k=1:nIons

    sites(1).res={'Q72' 'Q94' 'V76' 'M90' 'W255' 'M90' 'W255' 'V293' 'F417'};
    sites(2).res={'F417' 'S416' 'Q72'};
    sites(3).res={'S69' 'Q72' 'Y144'}; % we'll leave this Na2 tho empty.

    nSites=numel(sites(k).res);
    sites(k).dirs=zeros(1,nSites);

    pt=ptrsI(k);
    ionLoc=[pd.X(pt); pd.Y(pt); pd.Z(pt)];

    atomLocs=zeros(3,nSites);
    sites(k).ptrs=zeros(1,nSites);
    sites(k).dist=zeros(1,nSites);
    sites(k).ionLoc=ionLoc;

    disp(' ');
    disp(['Ligand: ' ligandLabels{k} ' coords: ' num2str(ionLoc')]);
    disp(' Site atoms, distances and rel coords:')
    for j=1:nSites % loop over coordinating residues
        siteLines=find(pd.resNum==str2double(sites(k).res{j}(2:end)));
        sL=siteLines; % pdb file lines of every atom in the residue
        nLines=numel(siteLines);
        nm=pd.atomName(sL);
        siteLocs=[pd.X(sL); pd.Y(sL); pd.Z(sL)];
        dists=sqrt(sum((siteLocs-repmat(ionLoc,1,nLines)).^2,1));
        [sites(k).dist(j),atomInd]=min(dists); % find the index of the closest atom
        sites(k).ptrs(j)=sL(atomInd); % pdb line of the closest atom
        atomName=nm{atomInd};
        atomLocs(:,j)=siteLocs(:,atomInd);
        disp([num2str(sites(k).ptrs(j)) ' ' sites(k).res{j} ' ' atomName ' ' ...
            num2str(sites(k).dist(j)) '   ' num2str(atomLocs(:,j)'-ionLoc')]);
    end;

end;

