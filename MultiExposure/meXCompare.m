% meCompareVesicleCoords

[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[setPath infoPath]=ParsePath(pa);
[dirPath localPath]=ParsePath(setPath);

if ~iscell(fname)
    fname={fname};
end;

cd(setPath);
nfiles=numel(fname);
mis=cell(nfiles,1);
for fileIndex=1:nfiles
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    mis{fileIndex}=mi;
end;
%  last file is the merged image
mi=mis{nfiles};
nfiles=nfiles-1;

%% Find the distances and outliers
maxR=2;  % outlier limit
refInd=4;

nv=numel(mis{nfiles}.vesicle.x);
ok=true(1,nv);
dr=zeros(2,nv,nfiles);
r0(1,:)=mi.vesicle.x;
r0(2,:)=mi.vesicle.y;
for i=1:nfiles
    dr(1,:,i)=mis{refInd}.vesicle.x-mis{i}.vesicle.x;
    dr(2,:,i)=mis{refInd}.vesicle.y-mis{i}.vesicle.y;
    ok=ok & sqrt(sum(dr(:,:,i).^2))<maxR;
end;
rs=sqrt(sum(dr(:,ok,:).^2));

% Plot the movements
vmag=60;
scaleBarLength=10;
colors=[ 1 0 .6; .8 .6 0; .3 .6 1 ; .5 .3 0];
subplot(1,1,1);
set(gca,'colororder',colors,'nextplot','replacechildren');

% hist(rs(:),80);
% create vector plot
    vxs=NaN(3*nv,nfiles);
    vys=NaN(3*nv,nfiles);
for k=1:nfiles
    for i=1:nv
        if ok(i)
            vxs(i*3-2,k)=r0(1,i);
            vxs(i*3-1,k)=r0(1,i)+dr(1,i,k)*vmag;
            vys(i*3-2,k)=r0(2,i);
            vys(i*3-1,k)=r0(2,i)+dr(2,i,k)*vmag;
        end;
    end;
end;
    plot(vxs,vys,'-',vxs(1:3:3*nv,1),vys(1:3:3*nv,1),'kd','markersize',2);
    hold on;

nx=mi.imageSize(1);
x0=nx*.04;
x1=x0+scaleBarLength*vmag/mi.pixA;
y0=x0;
plot([x0 x1],[y0 y0],'k-','linewidth',2);
text((x0+x1)/2,y0+nx*.002,[num2str(scaleBarLength) ' Å shift'],...
    'verticalalignment','bottom','horizontalalignment','center');
legString=cell(nfiles,1);
for i=1:nfiles
    legString{i}=[num2str(round(mis{i}.ctf.defocus*10)/10) 'um'];
%     legString{2*i}=' ';
end;
legend(legString);
hold off;
xlabel('x pixel');
ylabel('y pixel');
title([localPath infoPath fname{nfiles+1}],'interpreter','none');

