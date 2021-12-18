% rsRefineMembraneModel3
% Newer version with Relion compatibility.
% fs 7 Dec 2021
% Use linear least-squares fitting to deduce a membrane profile from the
% vesicles in an image.  We use the vesicles chosen by rsVesicleChooser
% by looking at mi.vesicle.ok(:,4) as the selected ones.
%
% Given an info structure mi and a merged image, use
% the stored vesicle centers and radii to create basis functions of shells
% 1 pixel thick, which are then used in the fitting.


ds=4;           % downsampling relative to raw image
forceNewFiles=1;
symmetrize=1;   % force a symmetrical profile
useAllVesicles=0;  % otherwise use the vesicles selected by rsVesicleChooser

% Modifications to the ctfs
bFactor=80;
% alpha0=.02;

writeMiFiles=0;

% Fitting
modelHalfwidth=70;   % angstroms
fitDC=1;        % also fit a constant offset to the image.

figure(1); clf; SetGrayscale;

% If we've already run this, we don't ask for a filename
if ~(exist('fname','var') && exist('rootPath','var')) || forceNewFiles  % not already defined
    [fname pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    %     [fname pa]=uigetfile('*mi.txt','Select mi file');
    % pa
    [rootPath infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
end;

cd(rootPath);
nFiles=numel(fname);
% nfiles=1;%%%%
%%

mis=cell(nFiles,1);
mi1s=cell(nFiles,1); % to receive the modifed mi files.

tic
for fileIndex=1:nFiles
    disp(['Reading ' fname{fileIndex}]);
    mi=ReadMiFile([infoPath fname{fileIndex}]);
    %%
    %     iname=[mi.procPath mi.baseFilename 'm.mrc'];
    %     disp(['Reading image: ' iname]);
    m0=meLoadNormalizedImage(mi,mi.padImageSize);
    ds0=1;
    n=mi.padImageSize/ds;
    m=Downsample(m0,n);
    subplot(2,2,1);
    imags(GaussFilt(m,256/n(1)));
    drawnow;

    if fileIndex==1
        imgs=m;
    else
        imgs(:,:,fileIndex)=m;
    end;


    if ~isfield(mi,'tempPath')
        mi.tempPath='Temp';
    end;
    if ~exist(mi.tempPath,'dir')
        mkdir(mi.tempPath);
    end

    if useAllVesicles
        mi.vesicleModel=Crop(ones(27,1),33);
        mi.vesicle.ok= mi.vesicle.r<150 & ~isnan(mi.vesicle.s(:,1));

        numOk=sum(mi.vesicle.ok(:,3))
        subplot(2,2,2);
        imags(GaussFilt(mres,256/n(1)));
        vesIndices=find(mi.vesicle.ok(:,3));
        mi.vesicle.ok(:,4)=mi.vesicle.ok(:,3); % copy the choices
    else % use only the selected ones
        vesIndices=find(mi.vesicle.ok(:,4));  % true means a picked vesicle.
        nv=numel(vesIndices);
        disp([num2str(nv) ' vesicles to be fitted']);
    end;

    mis{fileIndex}=mi;

    %

    % Create an mi structure for modification
    mi1=mi;
    %         Force all the vesicles to the same amplitude
    mi1.vesicle.s=ones(size(mi1.vesicle.s,1),1)*median(mi.vesicle.s(find(mi.vesicle.ok(:,1)),1,1),'omitnan');
    % Make basis images
    nvs=round(modelHalfwidth/(ds*mi.pixA));
    nxs=nvs*ds;  % halfwidth of the oversampled model.
    nm=2*nvs+1;  % number of points in the model
    nx=2*nxs+1;
    if symmetrize
        nt=nvs+1;  % number of terms in the model fit
    else
        nt=nm;
    end;

    %
        if fileIndex==1
            shells=zeros([n nt nFiles]);   % basis without CTF
    end;

    disp(['Constructing ' num2str(nt) ' shell images.']);
    for i=1:nt % For each term, make the vesicle model vm one (or two symmetric) delta functions.
        j=(i-1-nvs)*ds+nxs+1;
        vm1=zeros(nx,1);
        vm1(j)=1;  % scale it according to thickness
        if symmetrize
            vm1(nx-j+1)=1;
        end;
        mi1.vesicleModel=vm1*ds;  % removed the mi.pixA factor
        %         mi.vesicleModel=vm1*mi.pixA;
        %         r0=meMakeModelVesicles(mi,n*ovs,vesIndices,0,0); % no filtering
        r0=meMakeModelVesicles(mi1,n,vesIndices,0,0); % no filtering
        shells(:,:,i,fileIndex)=r0; % ith model image

        q=Crop(r0,256);
        subplot(2,2,4);
        imacs(q);
        title(i);
        subplot(2,2,3);
        plot(sect(q));
        drawnow;
    end;
end; % for fileIndex
disp('All the shell bases are made.');
toc

% We now have the assembled variables
%  mis, imgs, shells; plus nFiles. We assume pixA and ds are the same.
%%
%
% -----given the shell images, do the fitting------

alpha0=.015;
    bFactor=100;
n1d=prod(n);

disp('Applying CTF to shells');
Rtot=zeros(n1d, nFiles, nt);
for fileIndex=1:nFiles
    mi1=mis{fileIndex};
    for i=1:numel(mi1.ctf)
        mi1.ctf(i).alpha=alpha0;
        mi1.ctf(i).B=bFactor;
        mi.ctf(i).alpha=alpha0;
    end;

    H0=ifftshift(meGetEffectiveCTF(mi1,n,ds)); % This is slow, so we compute it once.
    subplot(2,2,2);
    plot((0:n(1)/2-1)/(n(1)*ds*mi.pixA),sectr(fftshift(H0)));
    xlabel('Spatial frequency');
    ylabel('CTF')
    drawnow;
    for i=1:nt
        r=real(ifftn(fftn(shells(:,:,i,fileIndex)).*H0));  % Filter with the ctf
        Rtot(:,fileIndex,i)=r(:);
    end;
end; % for fileIndex

Rtot=reshape(Rtot,n1d*nFiles,nt);
Mtot=reshape(imgs,n1d*nFiles,1);

disp('Doing the least-squares')
niter=1;
% Do the actual least-squares here
vAccum=zeros(n);
disp('Least squares')
for j=1:niter
    if fitDC
        dcNorm=mean(Rtot(:,1));
        Rtot(:,nt+1)=dcNorm*ones(n1d*nFiles,1);  % constant term
        a=LinLeastSquares(Rtot,Mtot);
        disp(['Const val is ' num2str(a(nt+1)*dcNorm)]);
        a(nt+1)=[];  % remove the constant
        Rtot(:,nt+1)=[];
    else
        a=LinLeastSquares(Rtot,Mtot);
    end;
    %     Get the fitted models
    aFit=reshape(Rtot*a,[n nFiles]);
end;
figure(2);
nr=ceil(sqrt(nFiles));
nc=ceil(nFiles/nr);
for i=1:nFiles
    mysubplot(nc,nr,i);
    imags(aFit(:,:,i))
end;
figure(3);
for i=1:nFiles
    mysubplot(nc,nr,i);
    imags(imgs(:,:,i)-aFit(:,:,i));
end;
%
disp('Fix the model')
if symmetrize
    ax=[a ; flipud(a(1:nt-1))];  % full profile
else
    ax=a;
end;
figure(1);
subplot(2,2,3);
plot(ax,'.-','markersize',10); hold on;
plot(0*ax,'k-'); hold off;
a0=a;
asym=a;
%
% upsample the model
amu=meDownsampleVesicleModel(ax,1/ds);
subplot(224);
plot(amu);
%
return

% Try subtracting the model
mi.vesicleModel=amu;
%         r1=meMakeModelVesicles(mi,n,vesIndices); %  filtering on
%
% subplot(2,2,4);
% msub2=mres-r1/2;
% imacs(GaussFilt(msub2,128/n(1)));
% drawnow;
%
%

% Try to force zero baseline.
pedpt=0;  % number of points remaining before zero crossing
a=a0;
pedestal=ones(nt,1);
spike=zeros(nt,1);
a1=a;
if ~symmetrize
    a1(ceil((nt+1)/2):nt)=1;  % Force the 2nd half to be >0
end;
p1=find(a1<=0,1,'last');
if numel(p1)<1 || p1<1+pedpt
    p1=1+pedpt;
end;
pedestal(1:p1-pedpt)=0;
spike(p1-pedpt)=1;
if symmetrize
    %         pedestal(np-p1+2:np)=0;
    %         spike(np-p1+2)=1;
else
    a2=a;
    a2(1:ceil((nt+1)/2))=1;  % Force first half ot be positive
    p2=find(a2<=0,1,'first');
    if numel(p2)<1 || p2>nt-1-pedpt
        p2=nt-1-pedpt;
    end;
    pedestal(p2+pedpt:nt)=0;
    spike(p2+pedpt)=1;
end;
a=a.*pedestal;
% R1=[R*a R*pedestal R*spike];
R1=[R*a R*pedestal];
a1=LinLeastSquares(R1,mres);

% d1=[a pedestal spike]*a1;
d1=[a pedestal]*a1;
if symmetrize
    d1=[d1; flipud(d1(1:nt-1))];
end;
d1x=meDownsampleVesicleModel(d1,1/ds);
d1x=circshift(d1x,-1);

ndx=numel(d1x);
ndxc=ndx-2*(ds*(p1-1)+1);
d1xc=Crop(d1x,ndxc);
d1xc=d1xc*2/max(d1xc);  % force maximum to be 2.
subplot(2,2,1);
plot(1:ndxc,d1xc,'.-',(ndxc+1)/2,d1xc((ndxc+1)/2),'r+','markersize',10);  % final composite profile
xlabel([num2str(mi.pixA) ' A per pixel']);
drawnow;

outName=[mi.tempPath mi.baseFilename 'ds' num2str(ds) 'sym' num2str(symmetrize)...
    'all' num2str(useAllVesicles) 'B' num2str(bFactor)...
    'alph' num2str(alpha0) 'nFiles' num2str(nFiles)];
print('-djpeg','-r300',[outName '.jpg']);  % save the plots as a jpeg.
disp('Figure image:');
disp(outName);

vm.vesicleModel=d1xc;
vm.pixA=mi.pixA;
save([outName '.mat'],'vm');
if writeMiFiles
    for i=1:nFiles
        mis{i}.vesicleModel=d1xc;
        theName=WriteMiFile(mis{i});
        disp(['Wrote ' theName]);
    end;

end;
