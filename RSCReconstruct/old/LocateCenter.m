% compute all the cc for different references and compute the p_img
% need to separate alphpa and the other angles
% inputs: templates(refs),rotImgs(cropped imgs),angleList(the other
% angles), alphas
% meaningful returns: pcoords(peak coordinates) and refi (index for ref and angle)

figure(4);
SetGrayscale;

ref=circshift(rsMakeTemplates([0 30 0],map),trans);
rotRefs=rsRotateImage(ref,alphas);

nref=size(angleList,1); % no. of refs
nimg= size(rotImgs,3); % no. of images
ns=size(rotImgs,1); % size of image 
nrs=size(templates,1); % size of ref, should be the same as ns for cc
nk=8; % crop out the central part for likelihood comparison
disp('Computing cross-correlations');
ccrs=single(zeros(ns,ns,nref));
ccim=single(zeros(nk,nk,nimg)); %cc coeff
ccri=single(zeros(nk,nk,nimg)); %cc reference index
unRotImgs=rotImgs;  % rotate the alpha back, may have better way
for i=1:nimg
    unRotImgs(:,:,i)=rsRotateImage(rotImgs(:,:,i),-alphas(i));
end;
%frefs=fftshift(fftn(templates),3);
frefs=zeros(ns,ns,nref);
for j=1:nref
    frefs(:,:,j)=fftn(templates(:,:,j));
end;
%fimgs=fftshift(fftn(rotImgs),3);  this is a wrong way to calculate fftn
%for a stack of image
%fimgs=fftshift(fftn(rotRefs),3);
%fimgs=fftn(rotImgs);
%fimgs=fftn(rotRefs);
pcoords=zeros(2,nimg);
refi=zeros(nimg,1);
refi2=refi;
coef=zeros(nref,1);
% cc is scaled by the norm
for i=1:nimg
    fms=fftn(unRotImgs(:,:,i));
    %fms=fimgs(:,:,i);
    for j=1:nref
        coef(j)=norm(templates(:,:,j));
        ccrs(:,:,j)=real(fftshift(ifftn(fms.*conj(frefs(:,:,j)))))/(coef(j)^2);
        %ccrs(:,:,j)=real(ifftn(fms.*conj(frefs(:,:,j))));
        %imacs(fftshift(ccrs(:,:,j)));
        %drawnow;
    end;
    [ccmx ccmi]=max(ccrs,[],3);  % ccmi is the index for angleList and templates
    ccim(:,:,i)=Crop(ccmx,nk); % ccmx is the for cc map for each image
    ccri(:,:,i)=Crop(ccmi,nk); % cc reference index
    [presentmax,jx,jy]=max2di(ccmx.*fuzzymask(n,2,n*0.05,0)); % here is the coordinates for the peak in cc map
    % sth wrong, jx jy is not at the center, apply fuzzymask didn't help
    refi(i)=ccmi(round(jx),round(jy)); % refs index
    
    % here compute the p_img
    
    %plotting
    subplot(221);
    %imacs(rotImgs(:,:,i));
    imacs(rotRefs(:,:,i));
    title(['Img ' num2str(i)]);
    subplot(222);
    imacs(ccim(:,:,i));
    pcoords(:,i)=[jx jy];
    title(['CC peak ' num2str([jx jy])]);
    subplot(223);
    imacs(ccmi);
    title('Reference index');
    %text(round(jx),round(jy),'Mid');
    %title(num2str(refi(i)));
    subplot(224);
    %imacs(templates(:,:,ccmi(24,24)));
    %imacs(templates(:,:,refi(i)));
    imacs(ccri(:,:,i));
    title(['Ref ' num2str(refi(i))]);  % center
    refi2(i)=ccmi(24,24);
    %title(num2str(refi(i)));
    drawnow;
end;

%% here compare the ref and image


figure(5);  % now draw the templates and the matche templates
SetGrayscale;
for i=1:nimg
    subplot(4,nimg,i);
    imacs(templates(:,:,refi2(i))); %
    axis off;
    title(['Ctr' ' ' num2str(refi2(i))]);
    subplot(4,nimg,i+nimg);    
    imacs(rotRefs(:,:,i));
    axis off;
    title('R');
    subplot(4,nimg,2*nimg+i);
    imacs(rotImgs(:,:,i));
    axis off;
    title('I');
    subplot(4,nimg,3*nimg+i);
    imacs(templates(:,:,refi(i)));
    title(['Max' ' ' num2str(refi(i))]);
end;


% %%  start to compute the likelihood p_img take sigmaN
% figure(6);
% SetGrayscale;
% % for one point, the cc peak
% % p_img=zeros(nimg,1);
% % for i=1:nimg
% %     ref=templates(:,:,refi(i));
% %     ref_t=circshift(ref,round(pcoords(:,i))'-[24 24]);
% %     %ref_t=circshift(ref,-1,1);  or compute a 2-D matrix of 2-D ref_t around the
% %     %image center/particle center, p_img will be a stack of 2-D matrix say 8*8*nimg
% %     p_index=(norm(unRotImgs(:,:,i)-ref_t))^2/(2*sigmaN^2);
% %     p_img(i)=exp(p_index);
% % end;
% 
% % rewrite the previous part for a stack of matrix
% %nk=8; % the central part to compute likelihood 
% Q_imgs=zeros(nk,nk,nimg); % the log of p_imgs 
% q_img=zeros(nk, nk);
% iref=zeros(2,nimg);
% for i=1:nimg
%     for k=1:nk*nk
%         ref=templates(:,:,ccri(k));
%         ref_t=circshift(ref,round(pcoords(:,i))'-ns/2);
%         residual=unRotImgs(:,:,i)-ref_t;
%         imacs(residual);
%         q_img(k)=norm(residual)^2;
%         
%         
%     end;
%     [q_max,ik,jk]=max2d(q_img);
%     iref(:,i)=[ik,jk];
%     Q_imgs(:,:,i)=q_img;
%     subplot(2,nimg,i);
%     imacs(q_img);
%     title('q');
%     subplot(2,nimg,nimg+i);
%     imacs(templates(:,:,ccri(ik,jk)));
%     title('ref');
% end;

% p_c ~ N()?peak at center, considering all within fuzzymask, now take p_c as a delta function?
%p_c=exp(norm(coords(:,k)-[ns/2 ns/2])^2/2sigmaK);
%% normalize it with latent variables
%ptomax=p_img(for different reference center)*p_c*p(w|k l)*p(k,l)
