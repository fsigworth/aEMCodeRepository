% rlCheckDefParticles.m
% run after rlCheckDefocus to find representative particles
%  Note that our alpha0 = -90 - rlnAnglePsi for iso particles,
%                alpha0 =  90 - rlnAnglePsi for rso particles
cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2/');
% mergeDir='Merged_z/';
iOffset=1223;
% micros=1:100;
% for i=micros
%     j=i+iOffset;
%     mi=si.mi{j};
%     disp([num2str([i mnParts(j) mCls(j,2)],4) '   ' mi.baseFilename]);
% end;

%%
% index=input('Index? ');
siIndices=zeros(0,1);
siGood=false(0,1);
for i=1:10  % index of micrographs
    index=i+32;
    j=index+iOffset;
    mi=si.mi{j};
    cParts=cParticles{j};
%     clsParts=cls(cParts);
    if i==1
        acParts=cParts;  % si particle indices
    else
        acParts=[acParts; cParts];
    end;
end;
sParts=pInds(acParts);
clsParts=cls(acParts);

sOffset=min(sParts);
sMax=max(sParts);
% stName=[siName(1:end-5) 'tack.mrcs'];
% suName=[siName(1:end-7) 'tustack.mrcs'];
% mp=ReadMRC(stName,sOffset,sMax-sOffset+1);
% mu=ReadMRC(suName,sOffset,sMax-sOffset+1);
% imgs=mp(:,:,sParts-sOffset+1);
% imgu=mu(:,:,sParts-sOffset+1);
% save sq81_imgs.mat imgs imgu
load sq81_imgs.mat
%%
figure(6);
% rimgs=rsRotateImage(imgs,si.alpha0(sParts));
% imagsar(BinImage(GaussFilt(rimgs,.02,1),4),.005);
imagsar(BinImage(rimgs,4),.005);
mk=struct;
mk.marker='gs';
mk.markerSize=52;
j=0;
for ind=find(clsParts==2)'
    j=j+1;
    if j==15
        mk.marker='gs'; % markers for after first micrograph
    end;
    imagsar('mark',ind,mk);
end;
return

  