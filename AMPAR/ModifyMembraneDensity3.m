% ModifyMembraneDensity3.m
% Compute the density map from pdb coordinates,
% and modify TM density of the Kv or AMPAR channel by subtracting the 1D
% membrand model.
% This code is used to create models for picking and simulations.

model='Kv';  % Alternatives: AMPAR, Slo2.2
% model='Slo2';
doWrite=1;

switch model
    case 'Kv'
        pdbName='~/Structures/kv1.2/Kv1.2Tetramer.pdb';
        miName='~/Structures/kv1.2/KvMembraneModel_mi.txt';
        mctr=132; % membrane center is at y=130.
        oldMapName='/Users/fred/Structures/kv1.2/KvMap.mat';
        outMapName={'~/Structures/kv1.2/KvMapMbnSub.mat' '~/Structures/kv1.2/KvMapMbnAdd.mat'};
        outMapName2={'~/aEMCodeRepository/AMPAR/KvMap' '~/aEMCodeRepository/AMPAR/KvMapMbnAdd'};
        
        rotMat=eye(3);
    case 'AMPAR'
        pdbName='~/Structures/AMPAR/3KG2.pdb';
        miName='~/Structures/kv1.2/KvMembraneModel_mi.txt';
        mi=ReadMiFile(miName);
        mctr=41; % membrane center is at z=41A.
        outPixA=2.9;
        oldMapName='~/aEMCodeRepository/AMPAR/3KG2mapsub2.9A.mat';
        outMapName='~/aEMCodeRepository/Data/3KG2mapsub58.mat';
        outMapName2='~/aEMCodeRepository/Data/3KG2mapsub29.mat';
        rotMat=EulerMatrixInverse([90 2 -90]*pi/180)*EulerMatrixInverse([-48 35 -10]*pi/180);
    case 'Slo2'
        pdbName='~/Structures/Slo2.2/5a6e.pdb';
        
        miName='~/Structures/kv1.2/KvMembraneModel_mi.txt';
        mi=ReadMiFile(miName);
        mctr=105; % membrane center is at z=28A.
        outMapName='~/Structures/Slo2/Slo2AMbnSub.mat';
        outMapName2='~/aEMCodeRepository/AMPAR/Slo2AMbnSub';
        rm=zeros(3,3,4);
        for i=1:4
            rm(:,:,i)=EulerMatrixInverse([90*(i-1) 180 0]*pi/180);
        end;
        rotMat(2,:,:)=rm(3,:,:);
        rotMat(3,:,:)=rm(2,:,:);
        rotMat(1,:,:)=rm(1,:,:);
end;
[coords, types]=ReadPDBAtoms(pdbName);
%%
rotCoords=[];
rotTypes=[];
cShift=[-98 -98 0];
cOffset=repmat(cShift',1,size(coords,2));
for i=1:size(rotMat,3)
    rotCoords=[rotCoords rotMat(:,:,i)*(coords+cOffset)];
    rotTypes=[rotTypes ; types];
end;

tic
[totalDens, protDens]=SolventAndProteinDensity(rotCoords, rotTypes);
toc

solDens=totalDens-protDens-totalDens(1); % subtract the background water dens.
%  protein will be at -4.95V.
figure(1);
ShowSections(protDens);
%%



% switch model
%     case 'Kv'
        n=size(totalDens,1);  % it's 183
% %         n1=96;  % size at 2 A / pixel
        
        %     md0=mi.vesicleModel;
        md0=[  0.000857023     0.090249     0.183498     0.264517     0.317358     0.336271      0.34617     0.382163 ...
            0.479292     0.657588     0.891925      1.14217       1.3681      1.53902      1.65281      1.71687 ...
            1.73874      1.72836       1.7028      1.68149      1.68375      1.71868      1.76457      1.78947 ...
            1.76148      1.66201      1.51257      1.34799       1.2032      1.10801      1.07677      1.11869 ...
            1.2427      1.44227      1.66433      1.84037      1.90211      1.81423      1.64055      1.47786 ...
            1.42284      1.53123      1.73572      1.92806            2      1.88308      1.62857      1.32753 ...
            1.07118     0.921539      0.85281     0.810007     0.737981     0.599094     0.408389     0.198421 ...
            0.00199102 ]';
        % pixA of the membrane model is 1.247
        %     md=real(ifft(Cropo(fft(md0),71)));
        ds=1/1.247;  % scale to 1 angstrom
        md=meDownsampleVesicleModel(md0,ds)*ds*1.247;
        md=.5*(md+flipud(md));  % symmetrize
        diffDens=totalDens-totalDens(1);
        modDens=diffDens;
        mbnDens=zeros(n,n,n,'single');
        mbnDisc=fuzzymask(n,2,0.37*n,.05             *n);
        %     Subtract the membrane density within the molecular volume.
        solDens0=4.88;
        mbnScale=.3;
        ctr=(numel(md)+1)/2;
        for iy=1:numel(md)  % loop over the membrane y values
            y=mctr+iy-ctr;
            modDens(:,y,:)=diffDens(:,y,:)+solDens(:,y,:)*md(iy)*mbnScale;
            mbnDens(:,y,:)=repmat(md(iy)*mbnScale,1,n).*mbnDisc*solDens0;
        end;
        totDens=modDens+mbnDens;
        mbnOffsetA=mctr-ceil((n+1)/2);
        %     modDens=shiftdim(modDens,2);
        ShowSections(GaussFilt(Crop(totDens,160),.1),[81 51 91]);
        drawnow;
%% 
switch model
    case 'Kv'
        %  Downsample to 2 A
        %     modDens2=GaussFilt(Downsample(Crop(modDens,216),108),.2)*2;
        modDens2=Downsample(Crop(modDens,216),108)*2;
        modDens2=circshift(modDens2,[0 4 0]);  % shift by 4 pixels upwards
        totDens2=Downsample(Crop(totDens,216),108)*2;
        totDens2=circshift(totDens2,[0 4 0]);  % shift by 4 pixels upwards
        mbnOffsetA=mbnOffsetA+8;
        if doWrite  % Write out volume, and also skewed volume
            maps=modDens2;
            maps(:,:,:,2)=totDens2;
            for i=2

                map=maps(:,:,:,i);
                save(outMapName{i},'map','mbnOffsetA','pixA');
                disp(['saved: ' outMapName{i}]);
                map=shiftdim(map,2);
                save([outMapName2{i} '.mat'],'map','mbnOffsetA','pixA');
                disp(['saved: ' outMapName2{i}]);
            end;
        end;
            %             map=SkewVolume(map,36:40,.5);
%             save([outMapName2 'skew.mat'],'map','mbnOffsetA','pixA');
        %     Display it
        figure(4);
        % ShowSections2(GaussFilt(modDens2-sim,.05));
        q=shiftdim(GaussFilt(modDens2,.05)+GaussFilt(4*randn(108,108,108),.1),2);
        q=GaussFilt(shiftdim(totDens2,2),.2);
        ShowSections2(Crop(q,112),[55 55 75],45);
        nctr=57;
        subplot(336); hold on;
        plot(nctr,nctr+mbnOffsetA/2,'w+');
        hold off;
    case 'Slo2'
        modDens2=Downsample(Crop(modDens,160),80)*2;
        pixA=2;
        if doWrite  % Write out volume, and also skewed volume
            map=modDens2;
            save(outMapName,'map','mbnOffsetA','pixA');
            disp(['saved: ' outMapName]);
            map=shiftdim(map,2);
            save([outMapName2 '.mat'],'map','mbnOffsetA','pixA');
            WriteMRC(map,pixA,[outMapName2 '.mrc']);
            disp(['saved: ' outMapName2 ' as .mat and .mrc']);
%             map=SkewVolume(map,36:40,.5);
%             save([outMapName2 'skew.mat'],'map','mbnOffsetA','pixA');
         end;
        %     Display it
        figure(4);
        q=shiftdim(modDens2,2);
        
        nctr=55;
        ShowSections2(GaussFilt(Crop(q,104),.2),[nctr nctr nctr-10],45);
        subplot(335); hold on;
        plot(nctr,nctr+mbnOffsetA/2,'w+');
        hold off;
       
     case 'AMPAR'
        n=size(totalDens,1);  % about 207
        nvm=numel(mi.vesicleModel);
        md0=DownsampleGeneral(mi.vesicleModel,ceil(nvm*mi.pixA),mi.pixA);
        nvm0=numel(md0);
        hnvm0=ceil(nvm0/2);
        md1=zeros(1,1,n);
        md1(mctr-hnvm0+1:mctr-hnvm0+nvm0)=md0;
        mdMap=repmat(md1,[n n 1]);
        n1=NextNiceNumber(n/outPixA);
        d0=totalDens(1);  % ice density level
        %     increase the lipid level to 1.5x the model
        modDens=(totalDens-totalDens(1))+1.5*solDens/d0.*mdMap;
        modDens2=GaussFilt(DownsampleGeneral(modDens,n1,1/outPixA),.2);
        modDens2=circshift(modDens2,[0 -1 0]);
        ctr1=ceil((n1+1)/2);
        mctr1=round(mctr/outPixA);
        %     ShowSections(GaussFilt(modDens+mdMap,.05));
        
        modDens29=Crop(modDens2,80);  % expand it.
        mDown=6;
        modDens29=circshift(modDens29,[-1 0 -mDown]);
        mbnOffsetA=-(ceil((n+1)/2)-mctr+mDown*outPixA)
        
        ShowSections(modDens29,[38 38 41+round(mbnOffsetA/outPixA)]);
        %     ShowSections(modDens29);
        
        modDens58=Downsample(modDens29,40);
end; % switch

% if doWrite && ~kVMode
%     sim=single(modDens58);
%     pixA=outPixA*2;
%     save(outMapName,'sim','mbnOffsetA','pixA');
%     sim=single(modDens29);
%     pixA=outPixA;
%     save(outMapName2,'sim','mbnOffsetA','pixA');
% end;

return
%%
% Compare to our previous reconstructed image

% modDens1=DownsampleGeneral(modDens,96,1/(2*mi.pixA));

s=load(oldMapName);
figure(1);
ShowSections2(shiftdim(Crop(s.sim,112),2),[],45);

%
% figure(1);
% ShowSections2(GaussFilt(SharpFilt(modDens1,.12)+1.5*randn(n1,n1,n1),.17));

return
%%
map1Name='~/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96j/mrc/i19av01.mrc';
map2Name='~/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96j/mrc/i19bv01.mrc';
expMap=ReadMRC(map1Name)+ReadMRC(map2Name);
figure(2);
ShowSections2(expMap,[],45);



%%
% We use solDens as a mask.
if kVMode
    
    n1=size(TotalDens,1);
    n0=2*ceil(n1/(2*pixA));  % size of final map, even to avoid possibe bugs
    
    pixA=2.9;  % angstroms per pixel
    npix=64;
    
    
    cd('/Users/fred/aEMCodeRepository/AMPAR')
    [m,pixA]=ReadEMFile('KvMap.mrc');
    
    m=shiftdim(m,2);
    
    m(:,:,65:100)=m(:,:,65:100)*.5;
    % m(:,:,1:45)=0;
    ShowSections2(m);
    
    mbnOffsetA=48;
    pixA=2;
    
    map=m;
    if doWrite
        save KvMap.mat map pixA mbnOffsetA
    end;
end;
