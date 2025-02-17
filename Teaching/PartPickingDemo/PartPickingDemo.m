% PartPickingDemo.m
ds=4;
ind=3; % line to take in star file.
bsd=80; % box size on downsampled image.
templateFilt=.05;
hpFilt=.001;

cd /Users/fred/EMWork/Yangyu/20190325/MotionCorr/job061_micrographs_ctf_Topaz
[~,d]=ReadStarFile('micrographs_ctf.star');
d=d{1};
%%
ind=1;

[pa,nm,ex]=fileparts(d.rlnMicrographName{ind});
[m1,s]=ReadMRC([nm ex]);
n=NextNiceNumber(size(m1),5,8);
m1c=Crop(m1,n,0,mean(m1(:)));

ctfPars=rlCTFParsFromStruct(d,ind);
c1=CTF(n,ctfPars);
pixA=ctfPars.pixA;

nd=n/ds;
mDs=Downsample(m1c,nd);
ctfParsDs=ctfPars;
ctfParsDs.pixA=ds*ctfPars.pixA;
% % ctfParsDs.alpha=0;  % force a zero mean
pixAds=ctfParsDs.pixA;

figure(1);
subplot(111);
imaga(imscale(mDs,256,1e-3));
drawnow;

%% Mask out the very dark part.

msk=GaussFilt(GaussFilt(mDs<6.6,.1)>.01,.05);
q=(1-msk).*mDs;
cmsk=1-msk;
me=sum(q(:))/sum(cmsk(:));
md=(1-msk).*mDs+msk*me-me;
md=GaussHP(md,hpFilt*pixAds);
imags(md)
drawnow;

%%
% Pick up the 3D map
[map,sm]=ReadMRC('emd_7009_flat_0822A.mrc');
npd=size(map,1)/ds;
mapd=Downsample(map,npd);
mape=Crop(mapd,bsd).*fuzzymask(bsd,3,npd/2-2,4);

figure(2);
ShowSections(GaussFilt(mape,.2),[41 41 48]);

%%
% mdy=squeeze(sum(mape,2)); % sum along y
% mdz=squeeze(sum(mape,3));
cnd=CTF(bsd,ctfParsDs);
% mdyf=real((ifftn(ifftshift(cnd).*fftn(mdy))));
% figure(3);
% subplot(111);
% imags(Crop(mdyf,n/ds));

%%            %% Create the list of angles for the templates
ppVals.nGamma=3;
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

% Make the templates

disp(['Making ' num2str(nAngles) ' templates']);

% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),mape);
filtTemplates=SharpFilt(allTemplates,templateFilt*pixA*ds,0,1);
filtTempVars=zeros(nAngles,1);
normTemplates=zeros(bsd,bsd,nAngles,'single');
for i=1:nAngles
    q=filtTemplates(:,:,i);
    filtTempVars(i)=q(:)'*q(:);
    normTemplates(:,:,i)=q/sqrt(filtTempVars(i));
end;



% imagsar(filtTemplates);

% Make "junk" templates
partRadius=.6;
partFuzz=.1;
half=zeros(bsd,bsd);
half(:,round(bsd/2):bsd)=1;
half=SharpFilt(half.*fuzzymask(bsd,2,npd*partRadius,partFuzz*bsd),templateFilt*pixA*ds);
norm=half(:)'*half(:);
half=half/sqrt(norm);
% imags(half);
angles=0:10:359;
nh=numel(angles);
halfTemplates=rsRotateImage(half,angles);
nt1=nAngles;
nt2=nAngles+nh;
% Add the half-moon templates to the bandpass-filtered templates
filtTemplates(:,:,nt1+1:nt2)=halfTemplates;
disc=SharpFilt(fuzzymask(bsd,2,npd*partRadius*.8,partFuzz*bsd),templateFilt*pixA*ds);
nt2=nt2+1;
filtTemplates(:,:,nt2)=disc;

% do ctf filtering and normalization
cTemplates=zeros(bsd,bsd,nt2,'single');
for i=1:nt2
    q=real(ifftn(ifftshift(cnd).*conj(fftn(filtTemplates(:,:,i)))));
    norm=q(:)'*q(:);
    cTemplates(:,:,i)=q/sqrt(norm);
end;
figure(3);
imagsar(cTemplates);

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

vStack=reshape(ccStack,prod(nd),nt2);
[mxVals,mxInds]=max(vStack');
% [mxVParticles,mxIParticles]=max(vStack(:,1:nt1)');
% [mxVNon,mxINon]=max(vStack(:,nt1+1:nt2)');
% pindImg=reshape(mxIParticles,nd);
% pccImg=reshape(mxVParticles,nd);
% nindImg=reshape(mxINon,nd);
% nccImg=reshape(mxVNon,nd);

%%




nonScale=1.2;
        blankRadius=[100 150]/(pixA*ds);
        th0=.4;
        th1=5;

blankMasks=cell(2,1);
blankMasks{1}=(1-fuzzymask(2*ceil(blankRadius(1)+3),2,blankRadius(1),2));
blankMasks{2}=(1-fuzzymask(2*ceil(blankRadius(2)+3),2,blankRadius(2),2));

% imaga(imscale(md+2*padTemplates(:,:,20),256,.001));
% hold on;
% plot(nd(1)/2+1,nd(2)/2+1,'gs','markersize',msz);
% pause(2);
% hold off;
figure(4);

b=0;

while b~='q'
    skip=0;
    switch b
        case 'u' % up the bad scaling
            nonScale=nonScale*1.05;
        case 'i' % down with bad scaling
            nonScale=nonScale/1.05;
        case 'j' % increase threshold
            th0=th0*1.05;
        case 'k' % keep more
            th0=th0/1.05;
        otherwise
            skip=0;
    end;
    
    if ~skip
        stackScales=ones(nt2,1);
        stackScales(nt1+1:nt2)=nonScale; % scale-up of non-particle cc values
        
        msz=40; % marker size
        
        ccMax=inf;
        
        
        
        vStackScaled=vStack;
        for i=1:nt2
            vStackScaled(:,i)=vStack(:,i)*stackScales(i);
        end;
        
        [mxVals,mxInds]=max(vStackScaled,[],2);
        
        ccImg=reshape(mxVals,nd);
        indImg=reshape(mxInds,nd);
        
        ccImgWk=ccImg;
        
        picks=zeros(0,4);
        imaga(imscale(GaussFilt(md,.4),256,.001));
        
        % drawnow;
        marker='-';
        hold on;
        while ccMax>th0
            [ccMax,i,j]=max2d(ccImgWk);
            goodBad=1+(indImg(i,j)>nt1);
%             ccImgWk=ccImgWk.*(1-fuzzymask(nd,2,blankRadius(1+isBad),2,[i,j]));
            ccImgWk=Mask(ccImgWk,[i j],blankMasks{goodBad});
            %     ccImgWk=ccImgWk.*(1-fuzzymask(nd,2,blankRadius,0,[i,j]));
            if ccMax<th1
                if goodBad==2 % a bad particle
                    marker='ro';
                    plot(i,j,marker,'markersize',4);
                else
                    marker='gs';
                    plot(i,j,marker,'markersize',msz);
                    picks(end+1,:)=[i j ccMax indImg(i,j)];
                end;
            end;
            
            %     imags(ccImgWk);
            %     title(ccMax);
            %     drawnow;
        end;
        txt=num2str([th0 nonScale size(picks,1)],4);
        text(1,1,txt,'verticalalignment','bottom','horizontalalignment','left',...
            'fontsize',16,'color',[1 1 0]);
        
        hold off;
        
        
        drawnow;
    end; % if skip
    [ix,iy,b]=Myginput(1);
end; % while
disp('done');


