% StackExtractor Given a set of info files, extract particles, rotate to
% the standard position, and write *st.mrc and *stu.mrc files (stack and
% unsubtracted stack) into the Stack/ directory.  Also store *stall.mrc the
% merged stack from all images, and *tsi.mat which contains the si structure.
% Contrast is reversed so protein is white.

boxSize=64;  % Size of boxes to be extracted from merged images.

padSize=NextNiceNumber(boxSize*1.5);
padMask=fuzzymask(padSize,2,padSize*.48,.04);
figure(1);
SetGrayscale;

[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);
if ~exist('Stack','dir');
    mkdir('Stack');
end;

totalNParts=0;
np=1e5;  % preliminary stack size
si=struct;
si.miIndex=     uint16(zeros(np,1));
si.miParticle=  uint16(zeros(np,1));
si.alpha0=      single(zeros(np,1));
si.yClick=      single(zeros(np,1));  % in units of si.pixA
si.rVesicle=    single(zeros(np,1));
si.mi=cell(0);
nfiles=numel(fname);
ctfs=single(zeros(boxSize,boxSize,nfiles));
pixA0=0;  % unassigned value

for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);  % Load the mi file
    %     check if there is something to do.
    if ~isfield(mi,'stackPath')
        mi.stackPath='';
    end;
    mi.stackPath='Stack/';
    np=size(mi.particle.picks,1);
    if np>0
        iname=[mi.procPath mi.baseFilename 'm.mrc'];
        if FileExists(iname)
            m=ReadEMFile(iname);
        else
            error(['Could not find the merged image ' iname]);
        end;
    end;
    % Get the vesicle model mVes
    vesOk=0;
    if isfield(mi,'tempPath') && exist([mi.tempPath mi.baseFilename 'v.mrc'],'file')
        vfile=[mi.tempPath mi.baseFilename 'v.mrc'];
        disp(vfile);
        mVes=ReadMRC(vfile);
        vesOk=(size(mVes,1)==size(m,1));
    end;
    if ~vesOk
        disp('Make Vesicles');
        mVes=meMakeModelVesicles(mi,size(m));
    end;
    msub=m-mVes;
    n=size(m);
    ds=mi.imageSize(1)/n(1);  % downsample factor
    pixA=mi.pixA*ds;
    if pixA0==0
        pixA0=pixA;
    end;
    if pixA0~=pixA
        disp(['Change in pixA values: ' num2str([pixA0 pixA]) '  ' fname{fileIndex}]);
    end;
    %     Inverse filter
    if numel(mi.noiseModelCode)>0
        Ti=meGetNoiseWhiteningFilter(mi,n);  % for micrograph
        Tic=meGetNoiseWhiteningFilter(mi,boxSize,ds);  % for particle ctf
        msubFilt=real(ifftn(fftn(msub).*ifftshift(Ti)));
        mFilt=real(ifftn(fftn(m).*ifftshift(Ti)));
    else
        disp('Pre-whitening filter hasn''t been set.')
        mfilt=m;
        Tic=1;
    end;
    
    % ctf, including the prewhitening filter effects
    freqs=RadiusNorm(boxSize)/pixA;
    [coeffs ctf]=meComputeMergeCoeffs(freqs,mi.ctf,mi.doses);
    ctfs(:,:,fileIndex)=ctf.*Tic;
    
    %%
    subplot(2,2,1);
    imacs(BinImage(mFilt,4));
    nv=numel(mi.vesicle.x);
    vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
    vRadii=mi.vesicle.r(:);
    stack=single(zeros(boxSize,boxSize,np));
    stackUnsub=single(zeros(boxSize,boxSize,np));
    sumImg=zeros(boxSize,boxSize);
    sumImgSub=zeros(boxSize,boxSize);
    
    j=0;
    for i=1:np  % scan the particle coordinates
        type=mi.particle.picks(i,3);
        if (type==16 || type==32)  % valid particle types
            j=j+1;  % counter
            pCoords=mi.particle.picks(i,1:2);
            ind=mi.particle.picks(i,4);
            if ind<1  % no vesicle marked
                dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                [minDist ind]=min(dists);
                if minDist > vRadii(ind) % outside of membrane, check distance from mbn.
                    distm=abs(dists-vRadii);
                    [minDist ind]=min(distm);
                end;
                mi.particle.picks(i,4)=ind;  % insert our estimated vesicle index.
            end;
            %             Compute the angle relative to the normal
            
            pVec=pCoords-vCoords(ind,:);
            %             type
            %             pCoords
            %             vCoords(ind,:)
            %             ind
            alpha=atan2(pVec(2),pVec(1))-pi/2;
            
            pimg0Sub=-SharpFilt(ExtractImage(msubFilt,round(pCoords/ds+1),padSize).*padMask,.45);
            pimgSub=Crop(grotate(pimg0Sub,-alpha),boxSize);
            sumImgSub=sumImgSub+pimgSub;
            stack(:,:,j)=pimgSub;
            
            pimg0=-SharpFilt(ExtractImage(mFilt,round(pCoords/ds+1),padSize).*padMask,.45);
            pimg=Crop(grotate(pimg0,-alpha),boxSize);
            sumImg=sumImg+pimg;
            stackUnsub(:,:,j)=pimg;
            
            subplot(2,2,1);
            imacs(pimgSub);
            title(mi.baseFilename,'interpreter','none');
            subplot(2,2,3)
            imacs(sumImgSub);
            title(j);
            subplot(2,2,2);
            imacs(pimg);
            subplot(2,2,4);
            imacs(sumImg);
            drawnow;
            
            jt=j+totalNParts;  % pointer to overall stack index
            si.miIndex(jt)=     fileIndex;
            si.miParticle(jt)=  i;
            si.alpha0(jt)=      alpha*180/pi;
            si.yClick(jt)=      hypot(pVec(1),pVec(2))/ds;  % in units of si.pixA
            si.rVesicle(jt)=    mi.vesicle.r(ind)/ds;
        end;
    end;
    nparts=j;
    stack=stack(:,:,1:nparts);  % truncate the stack
    stackUnsub=stackUnsub(:,:,1:nparts);
    
    totalStack(:,:,totalNParts+1:totalNParts+nparts)=stack;
    totalNParts=totalNParts+nparts
    save([infoPath fname{fileIndex}],'mi');
    disp([infoPath fname{fileIndex} ' written']);
    si.mi{fileIndex}=mi;  % store a copy of the micrograph info
    
    
    outname=[mi.stackPath mi.baseFilename 'st.mrc'];
    WriteMRC(stack,pixA,outname);
    outname=[mi.stackPath mi.baseFilename 'stu.mrc'];
    WriteMRC(stackUnsub,pixA,outname);
    
end;
%%
si.pixA=pixA;
np=totalNParts;
% truncate the si arrays.
si.miIndex=     si.miIndex(1:np,1);
si.miParticle=  si.miParticle(1:np,1);
si.alpha0=      si.alpha0(1:np,1);
si.yClick=      si.yClick(1:np,1);
si.rVesicle=    si.rVesicle(1:np,1);
si.ctfs = ctfs;

outname=[mi.stackPath mi.baseFilename 'stall.mrc'];
WriteMRC(totalStack,pixA,outname);

% save the stackInfo structure
save([mi.stackPath mi.baseFilename 'tsi.mat'],'si');
