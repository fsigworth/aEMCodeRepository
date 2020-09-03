% PartPickingDemoPreprocessor.m


% PartPickingDemo.m
ds=6;%% downsampling factor
bsA=300; % box size in angstroms
psA=150  % particle size in angstroms
% symmetry=3;
% refMapName='emd_7009_flat_0822A_304.mrc';
symmetry=4;
refMapName= 


templateFilt=.05; % inverse angstroms
templateFilt=.033; % inverse angstroms
hpFilt=.002;      % inverse angstroms


dataHalf=1;

cd ~/data/20200603/

starInPath='CtfFind/job005/';
% mrcInPath='/Volumes/UMass Seagate Expansion Drive/Yangyu_UMASS_Krios_20200603/A-g3/alignedMRC/';
mrcInPath='alignedMRC/';
% matOutPath='TemplatePicker3/';
matOutPath='~/mini_drobo/Yangyu/20200603/TemplatePicker3/';
matOutPath='~/mini_drobo/Yangyu/20200603/TemplatePicker4/';
CheckAndMakeDir(matOutPath,1);

% Read the star file
% cd /Users/fred/EMWork/Yangyu/20190325/MotionCorr/job061_micrographs_ctf_Topaz/

starName=[starInPath 'micrographs_ctf.star'];
disp(['Reading ' starName]);
[~,d]=ReadStarFile(starName);
d=d{1};
nl=numel(d.rlnMicrographName);

% Calculate the pixel size from the first entry
cPars=rlCTFParsFromStruct(d,1);
pixA=cPars.pixA;
% Get our box size in pixels
bsd=NextNiceNumber(round(bsA/(ds*pixA)),7,8);
disp(['Box size is ' num2str(bsd) ' pixels, ' num2str(bsd*pixA*ds) ' A.']);
% Pick up the 3D map

[map,sm]=ReadMRC(refMapName);
disp(['Map box size is ' num2str(sm.pixA*size(map,1)) ' A']);
% bsd=bsdA/(ds*sm.pixA)
% npd=size(map,1)/ds;
mag=sm.pixA/(pixA*ds);
mapd=DownsampleGeneral(map,bsd,mag); % Don't use an intermediate volume larger than 1GB
mapMask=fuzzymask(bsd,3,bsd/2-3,4); % Nearly the maximum circular region
mapMean=mapMask(:)'*mapd(:)/sum(mapMask(:));
mape=mapMask.*(mapd-mapMean);
figure(2);
ShowSections(GaussFilt(mape,.2),round(bsd*[.5 .5 .65])+1);

%%-------Create templates----------
% Create the list of angles for the templates
ppVals.nGamma=5;
ppVals.nAlpha=10;
ppVals.nBeta=4;
symmetry=3;

gammaStep=360/(symmetry*ppVals.nGamma);

%         hemiAngles run from alpha=[0..360) and beta=[0..90)
[hemiAngles, angleInds]=rsListHemisphereAngles(ppVals.nAlpha, ppVals.nBeta);
nHemiAngles=size(hemiAngles,1);

angleList=zeros(ppVals.nGamma,nHemiAngles,3,'single');
for j=1:ppVals.nGamma;
    gamma=(j-1)*gammaStep;
    for k=1:nHemiAngles
        angleList(j,k,:)=[hemiAngles(k,:) gamma];
    end;
end;
nAngles=numel(angleList)/3;

% ---------Make the templates

disp(['Making ' num2str(nAngles) ' templates']);

% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),mape);
filtTemplates=SharpFilt(allTemplates,templateFilt*pixA*ds,0,1);
filtTempVars=zeros(nAngles,1);
normTemplates=zeros(bsd,bsd,nAngles,'single');
for i=1:nAngles
    p=filtTemplates(:,:,i);
    filtTempVars(i)=p(:)'*p(:);
    normTemplates(:,:,i)=p/sqrt(filtTempVars(i));
end;

%  imagsar(filtTemplates);
%
% Make "junk" templates
partRadius=.6; % fractio of particle size
psd=psA/(pixA*ds);
partFuzz=.1;
half=zeros(bsd,bsd);
half(:,round(bsd/2):bsd)=1;
junkMask=fuzzymask(bsd,2,psd*partRadius,psd*partFuzz);
half=SharpFilt(half.*junkMask,templateFilt*pixA*ds);
halfNorm=half(:)'*half(:);
half=half/sqrt(halfNorm);
% imags(half);
angles=0:10:359;
nh=numel(angles);
halfTemplates=rsRotateImage(half,angles);
nt1=nAngles;
nt2=nAngles+nh;
% Add the half-moon templates to the bandpass-filtered templates
normTemplates(:,:,nt1+1:nt2)=halfTemplates;
% Add one more template, a full disc
disc=SharpFilt(junkMask,templateFilt*pixA*ds);
nt2=nt2+1;
normTemplates(:,:,nt2)=disc/sqrt(disc(:)'*disc(:));

% allVars=zeros(nt2,1);
% for i=1:nt2
%     q=normTemplates(:,:,i);
%     allVars(i)=q(:)'*q(:);
% end;
% plot(allVars);
% return

%% --------loop over micrographs------------
figure(1);
startInd=1;
endInd=nl;
switch dataHalf
    case 1
        endInd=floor(nl/2);
    case 2
        startInd=floor(nl/2)+1;
    otherwise
        % use the full range.
end;
disp(['Processing ' num2str(endInd-startInd+1) ' micrographs.']);
disp(['startIndex= ' num2str(startInd)]);
for ind=startInd:endInd
    [pa,baseName,ex]=fileparts(d.rlnMicrographName{ind});
    mrcName=[mrcInPath baseName ex];
    disp([num2str(ind,6) '  ' mrcName]);
    if ~exist(mrcName,'file')
        disp('  ...not found. Skipping.');
        continue;
    end;
    [m1,s]=ReadMRC(mrcName);
    m1sz=size(m1);
    n=NextNiceNumber(m1sz,7,ds);
    m1c=Crop(m1,n,0,mean(m1(:)));
    offsets=[n(1)-m1sz(1) n(2)-m1sz(2)]/2; % origin offsets in padded image.
    
    ctfPars=rlCTFParsFromStruct(d,ind);
    % c1=CTF(n,ctfPars);
    pixA=ctfPars.pixA;
    
    nd=n/ds;
    mDs=Downsample(m1c,nd);
    ctfParsDs=ctfPars;
    ctfParsDs.pixA=ds*ctfPars.pixA;
    % % ctfParsDs.alpha=0;  % force a zero mean
    pixAds=ctfParsDs.pixA;
    
%     figure(1);
%     subplot(111);
    % imaga(imscale(mDs,256,1e-3));
    % drawnow;
    
    %% Mask out the very dark part. Actually, mask anything beyond 6*sd
    me=mean(mDs(:));
    sd=std(mDs(:));
%     msk=GaussFilt(GaussFilt(mDs<6.6,.1)>.01,.05);
    msk=GaussFilt(GaussFilt(abs(mDs-me)>6*sd,.05)>.01,.05);
    p=(1-msk).*mDs;
    cmsk=1-msk;
    me=sum(p(:))/sum(cmsk(:));
    md=(1-msk).*mDs+msk*me-me;
    md=GaussHP(md,hpFilt*pixAds);
    imags(md)
    title([num2str(ind) '  ' mrcName],'interpreter', 'none');
    drawnow;
    
    %%
    
    %%
    % mdy=squeeze(sum(mape,2)); % sum along y
    % mdz=squeeze(sum(mape,3));
    cnd=CTF(bsd,ctfParsDs);
    % mdyf=real((ifftn(ifftshift(cnd).*fftn(mdy))));
    % figure(3);
    % subplot(111);
    % imags(Crop(mdyf,n/ds));
    
    
    % do ctf filtering and normalization
    cTemplates=zeros(bsd,bsd,nt2,'single');
    for i=1:nt2
        p=real(ifftn(ifftshift(cnd).*conj(fftn(normTemplates(:,:,i)))));
        norm=p(:)'*p(:);
        cTemplates(:,:,i)=p/sqrt(norm);
    end;
    if ind==1
        figure(3);
        imagsar(cTemplates);
    end;
    
    %% Make full-size templates
    stackScales=ones(nt2,1);
    stackScales(nt1+1:nt2)=4;
    padTemplates=Crop(cTemplates,nd,1);
    fTemplates=zeros([nd nt2],'single');
    ccStack=zeros([nd nt2],'single');
    fmd=fftn(md);
    for i=1:nt2
        fTemplates(:,:,i)=conj(fftn(ifftshift(padTemplates(:,:,i))));
        %     fTemplates(:,:,i)=stackScales(i)*conj(fftn(ifftshift(padTemplates(:,:,i))));
        ccStack(:,:,i)=real(ifftn(fmd.*fTemplates(:,:,i)));
    end;
    
    borders=[nt1 nt2];  % ends of separate ccs
    nGroups=numel(borders);
    goods=1; % which ones are good particles
    
    vStack=reshape(ccStack,prod(nd),nt2);
    % [mxVals,mxInds]=max(vStack');
    
    ccVecs=zeros(prod(nd),nGroups,'single');
    ccInds=zeros(prod(nd),nGroups,'uint16');
    
    border1=1;
    for i=1:nGroups
        
        [ccVecs(:,i),ccInds(:,i)]=max(vStack(:,border1:borders(i)),[],2);
        %     [mxVNon,mxINon]=max(vStack(:,nt1+1:nt2),[],2);
        border1=borders(i)+1;
    end;
    
    p=struct;
    p.ind=ind;
    p.baseName=baseName;
    p.micrographName=d.rlnMicrographName{ind};
    p.ds=ds;
    p.offsets=offsets;
    p.pixAds=pixAds;
    p.nd=nd;
    p.md=md; % image
    p.cnd=cnd; % ctf
    p.ctfParsDs=ctfParsDs;
    p.cTemplates=cTemplates;
    p.borders=borders;
    p.goods=goods;
    
    p.ccImgs=reshape(ccVecs,[nd nGroups]);
    p.ccInds=reshape(ccInds,[nd nGroups]);

    matName=[matOutPath baseName '.mat'];
    save(matName,'p');
    
end; % for ind

% % %%
% % 
% % 
% % 
% % 
% % nonScale=1.2;
% % blankRadius=[100 150]/(pixA*ds);
% % th0=.4;
% % th1=5;
% % 
% % blankMasks=cell(2,1);
% % blankMasks{1}=(1-fuzzymask(2*ceil(blankRadius(1)+3),2,blankRadius(1),2));
% % blankMasks{2}=(1-fuzzymask(2*ceil(blankRadius(2)+3),2,blankRadius(2),2));
% % 
% % % imaga(imscale(md+2*padTemplates(:,:,20),256,.001));
% % % hold on;
% % % plot(nd(1)/2+1,nd(2)/2+1,'gs','markersize',msz);
% % % pause(2);
% % % hold off;
% % figure(4);
% % 
% % b=0;
% % 
% % while b~='q'
% %     skip=0;
% %     switch b
% %         case 'u' % up the bad scaling
% %             nonScale=nonScale*1.05;
% %         case 'i' % down with bad scaling
% %             nonScale=nonScale/1.05;
% %         case 'j' % increase threshold
% %             th0=th0*1.05;
% %         case 'k' % keep more
% %             th0=th0/1.05;
% %         otherwise
% %             skip=0;
% %     end;
% %     
% %     if ~skip
% %         stackScales=ones(nt2,1);
% %         stackScales(nt1+1:nt2)=nonScale; % scale-up of non-particle cc values
% %         
% %         msz=40; % marker size
% %         
% %         ccMax=inf;
% %         
% %         
% %         
% %         vStackScaled=vStack;
% %         for i=1:nt2
% %             vStackScaled(:,i)=vStack(:,i)*stackScales(i);
% %         end;
% %         
% %         [mxVals,mxInds]=max(vStackScaled,[],2);
% %         
% %         ccImg=reshape(mxVals,nd);
% %         indImg=reshape(mxInds,nd);
% %         
% %         ccImgWk=ccImg;
% %         
% %         picks=zeros(0,4);
% %         imaga(imscale(GaussFilt(md,.4),256,.001));
% %         
% %         % drawnow;
% %         marker='-';
% %         hold on;
% %         while ccMax>th0
% %             [ccMax,i,j]=max2d(ccImgWk);
% %             goodBad=1+(indImg(i,j)>nt1);
% %             %             ccImgWk=ccImgWk.*(1-fuzzymask(nd,2,blankRadius(1+isBad),2,[i,j]));
% %             ccImgWk=Mask(ccImgWk,[i j],blankMasks{goodBad});
% %             %     ccImgWk=ccImgWk.*(1-fuzzymask(nd,2,blankRadius,0,[i,j]));
% %             if ccMax<th1
% %                 if goodBad==2 % a bad particle
% %                     marker='ro';
% %                     plot(i,j,marker,'markersize',4);
% %                 else
% %                     marker='gs';
% %                     plot(i,j,marker,'markersize',msz);
% %                     picks(end+1,:)=[i j ccMax indImg(i,j)];
% %                 end;
% %             end;
% %             
% %             %     imags(ccImgWk);
% %             %     title(ccMax);
% %             %     drawnow;
% %         end;
% %         txt=num2str([th0 nonScale size(picks,1)],4);
% %         text(1,1,txt,'verticalalignment','bottom','horizontalalignment','left',...
% %             'fontsize',16,'color',[1 1 0]);
% %         
% %         hold off;
% %         
% %         
% %         drawnow;
% %     end; % if skip
% %     [ix,iy,b]=Myginput(1);
% % end; % while
% % disp('done');
% % 
