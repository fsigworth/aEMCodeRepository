% TestTrimHole2.m
% Experiment with auto-trimming of thick carbon

mainPath='/Volumes/TetraData/EMWork/Hideki/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Merge/';
filename='DE_20120122_220714_719m.mrc';
mainPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Merged/';

mainPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Merged/';
infoPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Info/';
filename='004_sq02_1_04';
postscr='m.mrc';

% d=dir(mainPath);
% dirIndex=6;
% filename=d(dirIndex).name

figure(1);
SetGrayscale;

[m pixA]=ReadEMFile([mainPath filename postscr]);
load([infoPath filename 'mi.mat']);
n=size(m);
pars.n0=256;
pars.thresh=.6;
pars.width=.2;
pars.edge=.2;

emsk=AutoHoleMask(m,pixA,pars);

%% Save the mask in the mi structure

mask.merge='AND';
mask.encoding='RLE';  % RLE
mask.data=RunLengthEncode(emsk);

maskIndex=3;  % This is index for automasking
if ~isfield(mi,'mask')
    mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
    mi.mask(maskIndex)=mask;

save([infoPath filename 'mi.mat'],'mi');


%%
emsk=RunLengthDecode(mi.mask(3).data);
subplot(1,1,1);
emskx=ExpandImage(single(emsk),size(m,1)/(4*size(emsk,1)));
cMsk=repmat(rot90(emskx),[1 1 3]);
cMsk(:,:,3)=max(cMsk(:,:,3),.7);
cMsk(:,:,1:2)=max(cMsk(:,:,1:2),.3);
mscl=imscale(Downsample(m,n/4),256,1e-3)/256;
cImg=repmat(rot90(mscl),[1 1 3]);  % grayscale image
image(cMsk.*cImg);
