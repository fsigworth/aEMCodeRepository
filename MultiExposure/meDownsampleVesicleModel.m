function vd=meDownsampleVesicleModel(vm,ds)
% Interpolate the model vd0 and integrate over the vicinity of each output
% point, with the spacing between output points ds.
if ds==round(ds)  % an integer
    ovs=2;  % oversampling ratio for interpolation
else
    ovs=32;
end;
hwIn=floor(numel(vm)/2);
zds=zeros(ceil(ds*2),1);
vd1=[zds; vm(:); zds];  % pad with zeros to prevent getting out of range.
np1=numel(vd1);
hw1=floor(np1/2);  % half-width, i.e. points before center
% oversample
xsi=(1:1/ovs:np1)';  % should be (np1-1)*ovs+1 points.
vdi=interp1(vd1,xsi,'spline',0);

% Now integrate in the vicinity of each output point
hwo=ceil(hwIn/ds);  % output half-width
npo=2*hwo+1;
ctro=(hwo+1);
xso=(ctro-hwo*ds:ds:ctro+hwo*ds);  % original sample points
vd=zeros(hwo*2+1,1);  % output array
for i=1:npo
    xOrg=(i-ctro)*ds*ovs+hw1*ovs+1; % center is at hw1+1
    x0=round(max(xOrg-ovs*ds/2,1));
    x1=round(xOrg+ovs*ds/2);
    vd(i)=(sum(vdi(x0:x1))-(vdi(x0)+vdi(x1))/2)/(ovs*ds);
end;

% plot(-hwo:hwo,vd,'.-','markersize',10);
% axis([-20 20 0 2.5]);
% drawnow;
% 


% 
% 
% 
% if nargin<3
%     fc=.5;
% end;
% if ds>1  % downsample the model
%     nvp=numel(vd0);
%     fc1=fc/(min(4,ds));
%     vd1=GaussFilt([vd0;0*vd0],fc1);
%     hw1=floor(nvp/2);
%     hws=ceil(hw1/ds);  % halfwidth of downsampled copy
%     outPoints=(hw1+1-hws*ds:ds:hw1+1+hws*ds);
%     vd=interp1(vd1,outPoints','spline');
% 
% %     
% %     
% %     nv0=numel(vd0);
% %     vd0p=zeros(ceil(nv0*ds),1);  % a multiple of ds
% %     vd0p(1:nv0)=vd0;
% %     vds=circshift(vd0p,nv0+1-hw);  % Put the center in the middle
% %     vd=DownsampleGeneral(vds,2*hws+2,1/ds);
% %     vd(1)=[];  % make the number of elements odd
% %     vd=SharpFilt(vdf,.5/ds,.2/ds); % fourier filter
% %     vd=interp1(1:numel(vdf),vdf,(1+hw:ds:2*hws*ds+hw+1)','spline');
% %     vd=vdf(1+hw:ds:2*hws*ds+hw+1)*ds;
% elseif ds<1
% % upsample the model
%     nvp=numel(vd0);
%     hw1=floor(nvp/2);
%     hws=ceil(hw1/ds);  % halfwidth of upsampled copy
%     nv0=numel(vd0);
%     vd0p=zeros(nv0*2,1);  % pad by 2
%     vd0p(1:nv0)=vd0;
%     vds=circshift(vd0p,-hw1);  % ds*hws-hw would shift center to origin.
%     vdf=Downsample(vds,2*nv0/ds); % fourier filter
%     vdf=SharpFilt(vdf,.4,.2);
%     vd=circshift(vdf,hws);
%     vd=vd(1:2*hws+1);
% %     need to zero out the points that are zero on input...
% else
%     vd=vd0;
% end;
