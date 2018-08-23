% function [mi msub]=meRemoveVesicles(m, mi)
% Given a merged image m and the info structure mi, subtract vesicles from
% the merged image m according to the stored information in mi.  If this is
% called with no arguments, a file selector is put up to find the mi file.

% if nargin<1  % put up a file selector
[fname pa]=uigetfile('*mi.mat','Select an mi file');
cd(pa)
load(fname);
iname1=[mi.procPath mi.baseFilename 'mc.mrc'];
iname2=[mi.procPath mi.baseFilename 'm.mrc'];
if FileExists(iname1)
    m=ReadEMFile(iname1);
    iname=iname1;
elseif FileExists(iname2)
    m=ReadEMFile(iname2);
    iname=iname2;
else
    [iname pa]=uigetfile('*m.mrc','Find the merged image');
    m=ReadEMFile(iname);
end;

% Get image and pixel sizes
n=size(m,1);
ds=mi.imageSize(1)/n;  % downsampling factor of m
pixA=mi.pixA*ds;    % pixel size of m

%%
figure(1); clf; SetGrayscale;
model=zeros(n,n);
nv=numel(mi.vesicle.x);
model=meMakeModelVesicles(mi1,1:nv,n);
mv=m-model;
imacs(mv);

%%
[pai nmi ext]=fileparts(iname);

% Write the subtracted image back into the merged-image directory.
oname=[pai '/' nmi 'v.mrc'];
WriteMRC(mv,pixA,oname);
jname=[pai '/' nmi 'v.jpg'];
imwrite(uint8(imscale(mv)),jname);

