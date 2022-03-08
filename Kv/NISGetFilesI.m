function [names,pd,ptrsI,ligandLabels,sites,P]=NISGetFilesI;

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia/')
names.map='211021/cryosparc_P1_J340_005_volume_map_sharp.mrc';
names.pdb='211021/I_10-18-coot-21_real_space_refined_003 (1).pdb';
names.outDir='220216/';
names.dataName='NISMapDataI.mat'; % output file.

P=[35 0];  % rotation, z-translation

disp(['Loading the pdb file 5: ' names.pdb]);
pd=pdb2mat(names.pdb);

%% locate the ions
disp('Locating the ions in pdb 5');

NaPtrs=find(strcmp('Na',pd.element));
IPtrs=find(strcmp('IOD',pd.resName));
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

nIons=3;
sites=struct;
for k=1:nIons

    sites(1).res={'Q72' 'L413' 'F417' 'W255' 'V293' 'M90' 'Q94'};
    sites(2).res={'F417' 'S416' 'Q72'};
    sites(3).res={'S69' 'Q72' 'Y144'};

    % sites(k).dirs=[ 1     0       0      0     0      1      1];
    % siteAtom={'NE2' 'CD1'  'CZ'   'CZ2'  'CG2'  'CG'  'OE1'}; % all are closest
    %     except Q72:CG is closer than NE2 by .4A!
    nSites=numel(sites(k).res);
    sites(k).dirs=zeros(1,nSites);

    pt=ptrsI(k);
    ionLoc=[pd.X(pt); pd.Y(pt); pd.Z(pt)];

    atomLocs=zeros(3,nSites);
    sites(k).ptrs=zeros(1,nSites);
    sites(k).dist=zeros(1,nSites);

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

