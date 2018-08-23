% supposed we know the particle center in a cropped images, 

%%  start to compute the likelihood, take sigmaN
figure(9);
SetGrayscale;
ref=circshift(rsMakeTemplates([0 45 45],map),trans);   
rotRefs=rsRotateImage(ref,alphas);   % create non-noise image for debugging

nref=size(angleList,1); % no. of refs
nimg= size(rotImgs,3); % no. of images
ns=size(rotImgs,1); % size of image 
nrs=size(templates,1); % size of ref, should be the same as ns for cc
nk=8; % crop out the central part for likelihood comparison

for i=1:nimg
    %unRotImgs(:,:,i)=rsRotateImage(rotRefs(:,:,i),-alphas(i));
    unRotImgs(:,:,i)=rsRotateImage(rotImgs(:,:,i),-alphas(i));
end;
Q_imgs=zeros(nk,nk,nimg); % the log of p_imgs 
q_img=zeros(nk,nk);
iref=zeros(2,nimg);
ref_trans=zeros(2,nk,nk);
sqrmsk2=ones(nk,nk);

%sqrmsk2=-SquareWindow(nk,2);
%sqrmsk=-GaussMask(nk,2,4);
sqrmsk=1-Squaremask(nk,2,nk/4); % gauss index for clicking error

c_x=meshgrid(1:nk)-nk/2-1;
c_y=c_x';
ml=zeros(nref,1);
ccri=zeros(nk,nk);
for i=1:nimg
    unRotImg=unRotImgs(:,:,i);
    for k=1:nk*nk
        ref_trans(:,k)=[c_y(k),c_x(k)];
        %ml_threshold=0;
        %ml_ref=0;
        % compute all the refs
        subplot(3,nimg,2*nimg+i);
        for j=1:nref
            ref=templates(:,:,j);        
            ref_t=circshift(ref,ref_trans(:,k));
                 
            residual=unRotImg-ref_t;  % in this case don't need to scale the references
            
            
            %imacs(residual);
            %drawnow;

            ml(j)=-(norm(residual))^2/(2*sigmaN^2);
%             if ml>ml_threshold
%                 ml_threshold=ml;
%                 ml_ref=j;
%             end;
         end;
         [q_img(k),ccri(k)]=max(ml);
         
        %ref_trans(k)
    end;   
        
        %q_img(k)=ml_threshold;
        %ref(i)=ml_ref;

    
    [q_max,ik,jk]=max2d((q_img+sqrmsk).*sqrmsk2);
    ref_trans(:,ik,jk)
    [ik jk]
    %[q_max,ik,jk]=max2d(q_img);
   % iref(:,i)=-[ik,jk]+nk/2+1
    Q_imgs(:,:,i)=q_img;
    subplot(4,nimg,i);
    imacs((q_img+sqrmsk).*sqrmsk2);
    
    text(ik,jk,num2str(ccri(ik,jk)));
    text(3,5,num2str(ccri(3,5)));
    title(['q centered ' num2str(ccri(ref_trans(:,ik,jk)'+5)) num2str(ccri(3,5))]);
    subplot(4,nimg,nimg+i);
    imacs(templates(:,:,ccri(ik,jk)));
    title(['ref' num2str(ccri(ik,jk))]);
    subplot(4,nimg,i+2*nimg);
    residual=unRotImg-circshift(templates(:,:,ccri(ik,jk)),ref_trans(:,ik,jk));
    %residual=unRotImg-circshift(templates(:,:,ccri(ik,jk)),-iref(:,i));
    imacs(residual);
    title('residual');
    text(1,1,num2str(q_max));
    subplot(4,nimg,i+3*nimg);
    imacs(unRotImg);
    title('unrotImg');
    drawnow;
           
end;

% p_c ~ N()?peak at center, considering all within fuzzymask, now take p_c as a gauss function?
%p_c=exp(norm(coords(:,k)-[ns/2 ns/2])^2/2sigmaK);
%% normalize it with latent variables
%ptomax=p_img(for different reference center)*p_c*p(w|k l)*p(k,l)
