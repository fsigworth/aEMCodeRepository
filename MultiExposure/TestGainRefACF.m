% TestGainRefACF.m

ndis=64;

% Get two filenames
filter={'*.tif' '*.dm3'};
[fname fpath]=uigetfile(filter,'Select two blank-beam images','multiselect','on');
if isa(fname,'numeric') % user clicked cancel
    return
elseif isa(fname,'char')  % simple string
    numfiles=1;
else
    numfiles=numel(fname);  % cell array
end;
if numfiles<2
    error('Need two files for cross-spectrum');
end;
cd(fpath);
[pa nm ex]=fileparts(fname{1});
switch lower(ex)
    case '.dm3'
        m1=ReadDM3(fname{1});
        m2=ReadDM3(fname{2});
    case '.tif'
        m1=imread(fname{1});
        m2=imread(fname{2});
    otherwise
        error(['Invalid file extension ' ex]);
end;

mr1=RemoveOutliers(mr1);
mr2=RemoveOutliers(mr2);

corr=GetGainRefACF(mr1,mr2,6);
figure(1);
SetGrayscale;
subplot(2,2,1);
q=Crop(corr,ndis);
imacs(q);
subplot(2,2,2);
plot(q(:,ndis/2+1));

save RefCorrection corr
return


%%

ndis=64;
n=size(m1,1);
pre1='carbon_infocus';
m1=imread([pre1 '1.tif']);
m1=RemoveOutliers(m1);
me1=mean(m1(:));
m1=m1-me1;

m2=imread([pre1 '2.tif']);
m2=RemoveOutliers(m2);
me2=mean(m2(:));
m2=m2-me2;

ccraw=fftshift(real(ifftn(fftn(m1).*conj(fftn(m2)))))/(me1*me2);

bigcorr=Crop(corr,n);

subplot(2,2,3);
q0=Crop(ccraw,ndis);
q1=Crop(ccraw-bigcorr,ndis);
imacs(q0);
subplot(2,2,4);
plot([q0(:,ndis/2+1) q1(:,ndis/2+1)]);

% q1=q0-Crop(corr,ndis);
% imacs(q1);


