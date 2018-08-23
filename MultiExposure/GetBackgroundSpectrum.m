function mclfx=GetBackgroundSpectrum(mlf,mi0,ind)

procPath = '/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/10sep19aFavorites/Merge/';
baseName = '10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00002';
cd(procPath);
name=[baseName 'mlfMsk.tif'];
mm=imread(name);
% imacs(mm);
msk=single(mm==252);
msk=AffineTransform(msk,inv(mi0.mergeMatrix(:,:,ind)));
%
msk(:,1800:2048)=1;
msk(1950:2048,:)=1;
% imacs(msk);

bmsk=1-msk;
mclfx=mlf.*bmsk;

% make a self-consistent masking.
for i1=1:6
    bkf=GaussFilt(mclfx,.001);
%     imacs(bkf);
%     drawnow;
    mclfx=(mlf.*bmsk+msk.*bkf);
%     imacs(mclfx);
%     drawnow;
end;
% mclfx=mclfx-bkf;
