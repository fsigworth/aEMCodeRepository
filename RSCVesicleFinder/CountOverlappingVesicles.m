%CountOverlappingVesicles(mi)
% CheckVesicleOverlaps
nmi=100;
bins=0:.02:1;
minRadius=130;
maxRadius=400; % in pixels
crossHists=zeros(nmi,numel(bins));
fracs=zeros(nmi,1);
makeDiscs=1;
thk=2;

for imi=1:nmi
    mi=allMis{imi};
    mi1=mi;
    ds=8;  % downsampling from original image
    nv=numel(mi.vesicle.x);
    n=mi.padImageSize/ds;
    n2=prod(n);
    vesNorms=zeros(nv,1);
    vesImgs=zeros([n nv],'single');
    
    mi1.vesicle.ok(:,2)=mi1.vesicle.ok(:,2) & ...
        mi.vesicle.s(:,1)>0.5*median(mi.vesicle.s(:,1));
%         mi1.vesicle.r(:,1)>minRadius & mi1.vesicle.r(:,1)<maxRadius & ...
    goodVes=(mi1.vesicle.ok(:,2));
    nGood=sum(goodVes);
    goodVesInds=find(goodVes);
    
    mi1=mi;
    if makeDiscs
        v1=meMakeModelVesicleDiscs(mi1,n,goodVesInds,thk);
        fracs(imi)=sum(v1(:)>1)/sum(v1(:)>0.5);
        
        imags(v1);
        drawnow;
        
    else
        % Force a uniform amplitude for the vesicle cross-sections
        
        mi1.vesicleModel=2*fuzzymask(69,1,24,9);
        mi1.vesicle.s=0*mi1.vesicle.s;
        mi1.vesicle.s(:,1)=1;
        
        % disp('Making model vesicles');
        for i=1:nGood
            ind=goodVesInds(i);
            %     Compute cross-section images
            v1=meMakeModelVesicles(mi1,n,ind,0,0,1);
            vesImgs(:,:,i)=v1;
            vesNorms(i)=sqrt(v1(:)'*v1(:));
        end;
        % Get the cross-matrix
        % disp('Computing the cross-matrix');
        v2=reshape(vesImgs,n2,nv);
        crossNorms=vesNorms*vesNorms';
        cc=v2'*v2;
        cross=cc./crossNorms-eye(nv,nv);
        % imaga(cross*256)
        
        h=hist(cross(:),bins);
        fracs(imi)=sum(h(2:end))/sqrt(sum(h));
        crossHists(imi,:)=h;
        
    end;
    disp(num2str([imi,nGood, fracs(imi)]));
    
    % crossHist(1)=0; % sum over the whole histogram is nv^2
    % hist((mi.vesicle.s(find(goodVes),1)));
    % pause
end;


return


%% Here we create discs for each vesicle.


