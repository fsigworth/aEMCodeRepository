figure(1);
ax1=axes('position',[.1,.1,.8,.9]);
axes(ax1);
m=100*rand(256);
imshow(uint8(m),'parent',ax1);

iFH=imfreehand();
mask=createMask(iFH);

mc=repmat(m,1,1,3);
mc(:,:,1)=mc(:,:,1)+100*mask;

imshow(uint8(mc),'parent',ax1)

% pause(1);
% delete(iFH);  % -------workaround for deleting the imfreehand curve in 2015a

% -------Other workarounds that also do the job for me
% iFH2=findobj(ax1,'tag','imfreehand');
% delete(iFH2);
% % 
% % if numel(iFH2)>0
% %  set(iFH,'visible','off','hittest','off');
% % end;