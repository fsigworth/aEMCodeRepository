% rsCheckIntensity.m
% Using an allMis file, insert peak intensities into mi.ok(11:16)
% and median intensity of vesicle interiors in mi.ok(10), with 
% the fractional area coverage in mi.ok(9).

% The peak intensities are computed as the mean in the top percentiles x/3 to x,
% with x=10^y, y=[-3 -2.5 -2 -1.5 -1]

doWrite=1; % write out the modified allMis cell array.
% jpegDir='Picker_jpegs/';

% paths for working in /Volumes/D257/Hideki/201228/05_1_1/
inMisName='Picking_9/allMis_holemasked_i+f_7505.mat';
outMisName='Picking_9/allMis_holemasked_i2+f.mat';
imgPath='Merged_sms/';

checkVesicleMedian=1; % Do this calculation too.
minVesRadiusA=130;
maxVesRadiusA=400;


outMis=allMis;

ds1=8; % downsampled from 4x downsamlpled small image
ds2=16;
pcExponents=[-3 -2.5 -2 -1.5 -1];
npc=numel(pcExponents);
boostVal=30;

mMul=200; % image display scale
mMed=8;   % approx. image median value

figure(1);
set(1,'color',[.4 .4 .4]); % Gray background, no tools
set(1,'menu','none');
set(1,'toolbar','none');
subplot(224);
cla;

% inds=find(rej);
% inds=1:numel(fracs);
inds=1:numel(allMis);
nin=numel(inds);
disp([num2str(nin) ' indices']);
iptr=1;
b='0';
    intens=zeros(nin,npc);
    defoci=zeros(nin,1);
    nParts=zeros(nin,1);
    vesMedians=zeros(nin,1);
    %%
while iptr<nin
    
    ind=inds(iptr);
    mi=allMis{ind};
    defoci(ind)=mi.ctf.defocus;
    if isfield(mi,'particle') && isfield(mi.particle,'picks')
        nParts(ind)=size(mi.particle.picks,1);
    end;

    pcs=zeros(2,npc);
    spcs=pcs;
    micName=[imgPath mi.baseFilename 'ms.mrc'];
    if exist(micName,'file')
        m0=ReadMRC(micName)/mi.imageNormScale+mi.imageMedian;
        msk=meGetMask(mi,size(m0));
        fracUnmasked=sum(msk(:))/numel(msk);

        m0=m0.*msk+mi.imageMedian*(1-msk);
        n0=mi.padImageSize;
        n1=n0/ds1;
        n2=n0/ds2;
        %         m0c=Crop(m0,mi.padImageSize,0,mean(m0(:)));
        m1=GaussFilt(Downsample(m0,n1),.1);
        m2=GaussFilt(Downsample(m1,n2),.1);
        subplot(221);
        imaga((m1-mi.imageMedian)*mMul+128);
        title(ind);
        mx=m2(:);
        k=1;
        subplot(222);
        [mSort,iSort]=sort(mx,'descend');
        nmx=numel(mx);
        for j=1:npc
            if fracUnmasked<.1
                break
            end;
            pct=(10^pcExponents(j))*fracUnmasked;
            range=round(nmx*pct/3.1):round(nmx*pct);
            pcs(k,j)=mean(mSort(range));
            
            spcs(k,j)=std(mSort(range));
            mxMarked=mx;
            mxMarked(iSort(range))=mxMarked(iSort(range))+boostVal;
            
            mx2=reshape(mxMarked,sqrt(nmx)*[1 1]);
            if j==4 && k==1
                imags(mx2);
                title(num2str([k pct*100 pcs(k,j)]));
            end;
        end;
        intens(ind,:)=pcs(1,:);
        mi.ok(11:10+npc)=pcs(1,:)';
        
        if checkVesicleMedian && numel(mi.vesicle.x)>0
            % rsCheck(mi)
            
            minRadius=minVesRadiusA/mi.pixA;
            maxRadius=maxVesRadiusA/mi.pixA; % in pixels
            mi1=mi;
            mi1.vesicle.ok(:,2)=mi.vesicle.ok(:,2) & ...
                mi.vesicle.s(:,1)>0.5*median(mi.vesicle.s(:,1));
            goodVes=(mi1.vesicle.ok(:,2)) & mi1.vesicle.r(:,1)>minRadius ...
                & mi1.vesicle.r(:,1)<maxRadius;
            nGood=sum(goodVes);
            if nGood>0
                goodVesInds=find(goodVes);
                thk=50/(mi1.pixA*ds1);
                ns1.M=meGetImageScaling(mi.imageSize,n1,ds1);
                ns1.n=n1;
                v1=meMakeModelVesicleDiscs(mi1,ns1,goodVesInds,thk);
                vMsk=v1>.5 & meGetMask(mi1,n1);
                pMsk=vMsk(:);
                fracUnmasked=sum(pMsk)/numel(pMsk);
                m1x=m1(:);
                vesMedIntensity=median(m1x(pMsk));
                mi.ok(9)=fracUnmasked;
                mi.ok(10)=vesMedIntensity;
                m1pad=(m1-mi.imageMedian).*vMsk;
                subplot(223);
                imaga(m1pad.*mMul+128);
                title(vesMedIntensity);
                subplot(221);
                title([ind vesMedIntensity]);
                mi.ok(9:10)=[fracUnmasked vesMedIntensity];
                vesMedians(ind)=vesMedIntensity;
            end;
        end;
        
        outMis{ind}=mi;
    end; % if exist micName
    
        subplot(224);
        startInd=max(2,ind-5);
        plot(vesMedians(1:startInd-1),intens(1:startInd-1,4),'r.');
        hold on;
        plot(vesMedians(startInd:ind),intens(startInd:ind,4),'k.','markersize',10);
        hold off;
        axis(mMed+[-1 1 -1 1]);
        drawnow;
        
        disp(mi.ok([10 15 9])');
        
    iptr=iptr+1;
    %     else
    %         cla;
end;
%%
if doWrite
oldAllMis=allMis;
%%
allMis=outMis;
disp(['Writing ' outMisName '...']);
save(outMisName,'allMis');
disp('done.');
allMis=oldAllMis;
end;

return
%%
print('-djpeg','Intensity 210306.jpg');


%% Misc code for interactive viewing of results

% %         allMis{ind}.active=true;
% %     
% %         jpegName=[jpegDir mi.baseFilename '_i0.jpg'];
% %         %         tstr=sprintf('%04d file: %05d   def: %5.2f   res: %6.2f   FOM: %6.3f',...
% %         %             iptr,ind,mi.ctf.defocus, mi.ctf.resLimit,mi.ctf.ccc);
% %         if exist(jpegName,'file')
% %             im=imread(jpegName);
% %             subplot(121);
% %             image(im);
% %         else
% %     %         disp([num2str(ind) ' not found: ' jpegName]);
% %     %         imags(zeros(768));
% %             tstr='';
% %         end;
% % 
% % 
% % tstr=sprintf('%04d %4d %4d',iptr,round(pcs(1,1)*100), round(spcs(1,1)*100));
% % allMis{ind}.ok(11:numel(pcs)+10)=pcs(:);
% % subplot(121);
% % title(tstr,'fontsize',16,'fontweight','normal','color','w');
% % drawnow;
% % disp(tstr);
% % 
% % if b~='R'
% %     [x,y,b]=ginput(1);
% % end;
% % switch char(b)
% %     case 'i'
% %         tptr=MyInput('Pointer value ',iptr);
% %         iptr=tptr-1;
% %     case 'j' % junk
% %         allMis{ind}.active=false;
% %         fprintf(' X\n');
% %     case 'p' % go to previous
% %         iptr=iptr-2;
% %     case 'q' % quit (and set active=1...)
% %         break;
% % end;
% % iptr=iptr+1;
% % if iptr<1
% %     beep
% %     ind=1;
% % end;
% % % end;
% % 
% % disp('Done.');
% % return
% % %%
% % 
% % % Set the active flags
% % for i=1:nmi
% %     mi=allMis{i};
% %     mi.active=~rej(i);
% %     allMis{i}=mi;
% % end;
% % 
% % % save the allMisSel.mat
% % disp('Saving allMis8_sel.mat');
% % save allMis8_sel.m allMis
% % 
% % %%
% % % save the individual mi files
% % infoDir='Info8_sel/';
% % CheckAndMakeDir(infoDir,1);
% % for i=1:nmi
% %     mi=allMis{i};
% %     miName=[infoDir mi.baseFilename 'mi.txt'];
% %     WriteMiFile(mi,miName);
% %     if mod(i,1000)==0
% %         disp([num2str(i) ' ' miName]);
% %     end;
% % end;
% % 
% % % store all the mi files
