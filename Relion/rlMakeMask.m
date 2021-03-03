% rlMavkeMask
% Make the "sphere+ellipsoid" mask for Kv refinement
maskMode='sphere+ellipsoid';
% maskMode='sphere';
n0=112;
n1=size(maps,1);
nim=size(maps,4);
ctr=floor(n0/2+1);
switch maskMode
    case 'sphere+ellipsoid'
        msk0=min(1,fuzzymask(n0,3,[.32 .32 .15]*n0,.1*n0,[ctr ctr ctr+.18*n0])...
            +fuzzymask(n0,3,0.24*n0,0.12*n0,[ctr ctr ctr-.06*n0]));
        outName='MaskSphere+Elipsoid.mrc';
    case 'sphere'
        msk0=fuzzymask(n0,3,0.2*n0,0.08*n0,[ctr ctr ctr-.06*n0]);
        outName='MaskSphere.mrc';
end;

msk1=Crop(msk0,n1);
mskn=repmat(msk1,[1 1 1 nim]);
figure(1);
ShowSections((maps+.1).*repmat(msk1,[1 1 1 nim]),[],45)
figure(3);
ShowSections(msk1,[],45);
%         ShowSections(msk1,[],45);

% WriteMRC(msk1,pixA,outName);


return

% Make looser Kv mask

msk0=ReadMRC('KvMask2.26p128.mrc');
maps=ReadMRC('KvRef128.mrc');
msk1=GaussFilt(msk0>.1,.07);
ShowSections(msk1.*(GaussFilt(maps,.1)+10),[65 65 80],45);
WriteMRC(max(0,min(1,msk1)),2.26,'KvMask2.26wide.mrc');