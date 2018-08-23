% % GriddingReconTest
% Demonstration of using gridding functions to make projections at random
% angles and then do a 3D reconstruction from the projections.

n=48;  % Original volume size
ks=3;  % Gridding kernel size

% Create a real-space volume: two smooth spheres
rvol=single(fuzzymask(n,3,n*0.05,n*0.02,[n*.4 n/2+1 n/2+1]));
rvol=rvol+fuzzymask(n,3,n*0.2,n*0.02,[n/2+1 n*.7 n/2+1]);
figure(1);
ShowSections(rvol);
drawnow;

disp('3D FFT');
vol=gridMakePaddedFT(rvol);

disp('start projecting');
nim=50;  % number of random projections


comp=gridMakePreComp(n,ks);

refBetas=10:10:170;
refGammas=0:15:165;
refBetas=10:20:170;
refGammas=0:30:165;
nRefBetas=numel(refBetas);
nRefGammas=numel(refGammas);
k=0;
angs=[];
for beta=refBetas
    for gamma=refGammas;  % There are 6 gammas
        k=k+1;
        angs(k,:)=[0 beta gamma];
    end;
end;

angles=rsDegToEuler(angs);
nim=size(angles,1)

% nim=50;
% angles=zeros(nim,3);
% for i=1:nim
%     angles(i,:)=[rand rand*pi rand*2*pi];    
% end;
% 
% 


slices=complex(zeros(n,n,nim));
for i=1:nim
    p2=gridExtractPlane(vol,angles(i,:),ks);
    slices(:,:,i)=gridRecoverRealImage(p2,comp);
end;
figure(4);
% ImagicDisplay(slices);


disp('start reconstruction');

vol0=gridMakeNullFT(n,3);
vol1=vol0;
fslice0=gridMakeNullFT(n,2);

% We do Fourier insertions into two volumes.  vol0 just has constant planes
% inserted, and is used for normalization.  vol1 is the actual
% reconstruction volume.
% np=fslice0.np;
% np1=fslice0.np1;
fnorm=fuzzymask(n,2,n/2-3,2);
rnorm=fftshift(ifft2(ifftshift(fnorm)));
fslice0=gridMakePaddedFT(rnorm);
% fslice0.PadFT=fuzzymask(np1,2,np/2-3,2); % the constant plane is actually a fuzzy disc.
disp('Normalization reconstr');
tic
for i=1:nim
    vol0=gridInsertPlane(fslice0,vol0,angles(i,:));
end

figure(2);
ShowSections(abs(vol0.PadFT).^.3);  % Show the Fourier volume
drawnow;

disp('Actual reconstr');
for i=1:nim
    fslice=gridMakePaddedFT(slices(:,:,i));
    vol1=gridInsertPlane(fslice,vol1,angles(i,:));
end;
toc;
%% 

% Wiener filter "normalization"
epsi=1;

vol.PadFT=vol1.PadFT.*vol0.PadFT./(epsi+vol0.PadFT.^2);

% Get the reconstructed volume.
revol=gridRecoverRealImage(vol,comp);
figure(3);
ShowSections(revol);
