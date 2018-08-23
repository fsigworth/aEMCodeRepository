
nHalfTrans=4;
nt0=2*nHalfTrans+1;

dirs=[0 -1; -1 0; 0 1; 1 0]; % down, left, up, right
rings=zeros(nHalfTrans+1,8*nHalfTrans);

% hmRings=cell(nHalfTrans,1); % Treating the translation array as polar

pt=[1 1]*(nHalfTrans+1);  % central point
rings(1,1)=pt(1)+nt0*(pt(2)-1);
for i=1:nHalfTrans           % 'radius' away from central point
    side=2*i;                % length of side of square 'circumference'
    pt=pt+[1 1];  % move out by 1, to upper-right corner
    j0=1;
    for j=1:4 % directions
        for k=1:side % individual steps
            rings(i+1,j0)=pt(1)+nt0*(pt(2)-1);
            pt=pt+dirs(j,:);
            j0=j0+1;
        end;
    end;
end;
% rings=zeros(nHalfTrans+1,2*nHalfTrans);  % create an array with rows being rings
% for i=1:numel(hmRings);
%     nr=max(1,8*(i-1));
%     rings(i,1:nr)=(hmRings{i}(:,1)+nt0*(hmRings{i}(:,2)-1));
% end;




nRings=numel(hmRings);
cStep=floor(64/nRings);
colors=jet;
for ir=1:nRings;
    color=colors(ir*cStep,:);
    plot(hmRings{ir}(:,1), hmRings{ir}(:,2),'o','color',color);
    hold on
end;
hold off;

inImg=zeros(nt0,nt0);
ctr=(nt0+1)/2;
inImg(ctr,ctr:ctr+nHalfTrans)=rand(1,nHalfTrans+1);
inVec=inImg(:);

% same angle definition as rsRotateImage
for alpha=0:1:360

% outVec=reTransRotate(inVec,alpha);
nt2=numel(inVec);
nt0=sqrt(nt2);
nHalf=((nt0-1)/2);
rSizes=(1:nHalf)'*8;
shifts=round(rSizes*alpha/360);
outVec=inVec;
for i=1:nHalf
    ring=rings(i+1,1:rSizes(i));
    outVec(circshift(ring,shifts(i),2))=inVec(ring);
end;

outImg=reshape(outVec,nt0,nt0);
imacs(outImg);
drawnow;
end;