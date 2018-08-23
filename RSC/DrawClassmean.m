function DrawClassmean(classmeanI,angleList,lir,nos)
% classmeanI is m*m*n angleList(0,beta,gamma) lir(nimg*nalpha,nref)
ngamma=9;
%nref=size(angleList,1);
figure(13);
SetGrayscale;

if nargin<4
    nref=size(classmeanI,3);
    if nargin<3
        lir=zeros(nref,1);
    end;
    for i=1:nref    
        subplot(9,19,i);
        imacs(classmeanI(:,:,i));
        text(4,4,num2str(angleList(i,2)),'color','w');
        axis off;
        j=ceil(i/ngamma);
        if j<1
            j=1;
        end;
        title(num2str(round(max(lir(:,i)))));
        %title(num2str(round(lir(j))));
        %title(num2str(i));
        drawnow;
    end;
else
    nref=numel(nos);
    for  k=1:nref
        subplot(9,19,k);
        imacs(classmeanI(:,:,nos(k)));
        text(4,4,num2str(angleList(i,2)));
        drawnow;
    end;
end;
end