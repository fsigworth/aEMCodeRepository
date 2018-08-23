% PzxIterateMACIV2.m
%
% for one image I
%figure(9);
SetGrayscale;
subplot(111);
nimg= size(rotImgs,3); % no. of images
ns=size(rotImgs,1); % size of image 
%% build the RSC model and initiate image information
load 3KG2RotMap2.9A.mat
map2=Downsample(map,40)/2000;

map=Crop(map2,48);  % approx normalization to give s/n=10, pad to 48 pixels.
pixA=5.8;
fc=1/(40/pixA);
vol=SharpFilt(map,fc,fc/10);
ri.symmetry=2;
ri.angleStep=10;
ri.nGamma=round(90/ri.angleStep);
%ri.maxTranslation=15;
angleList=rsMakeTemplateAngles(ri);
%r=80*ones(nimg,1);
r=(1:nimg)'*2+50;   % radii for the image stack
dr=3; % The pixel offset of particle center from transmembrane region 

%betas=pi/2*ones(9,1)';
dbeta=pi/180*[-20:5:20,-20:5:0,0:5:15,20:-5:-20,3:20];  %  beta errors
%dbeta=0;
%betas=pi./[3,3,20000,3,3,20000,3,3,2,3,3,20000,3,3,20000,3,3,20000,3,3,20000,3,3,20000,3,3,20000];
%betas=pi./[6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2,6/4,3,2];
betas=pi./[6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2,6,3,2,2,2];  % true betas for a set of images
betas(1:27)=pi-betas(1:27); % inverse insertion

% % betas for another set
% betas=pi./[6,3,3,6,3,3,6,3,3,6,3,3,6,3,3,6,3,3,6,3,3,6,3,3,6,3,3]; 
betaMs=betas+dbeta;  % measured betas
betaMs=[betaMs pi*angleList(1:135,2)'/180];
%ralphas=reshape(-[alphas+5,60*ones(numel(alphas),1)',alphas-8],9,3)';
ralphas=zeros(nimg,1);
ralphas(1:45)=reshape(-[alphas+5,alphas,alphas-8,alphas-10,alphas+10],9,5)';  % measured -alphas with errors
ralphas(1:27)=pi+ralphas(1:27);
%ralphas=
y=(r+dr).*sin(betaMs');
y(1:27)=(r(1:27)-dr).*sin(betaMs(1:27)'); % measured y
%y=r*sin(60*pi/180);%.*sin(alphas); %set beta=60, alpha=0
nref=size(angleList,1);
SigmaN=1.5;
SigmaB=3;
SigmaR=1;
SigmaC=1;

lambdaE=EWavelength(200);
kWiener=.03;  % Wiener constant
symmetry=2;
lambda=0.5; % probability of right orientation
alphaRs=(-10:5:10)/180*pi;  % alpha map for pw, pwc
betaRs=(0:10:180)/180*pi; % beta map for pw, pwc

nalpha=numel(alphaRs);
nbeta=numel(betaRs);
ngamma=9;   % gamma number

d_a=10/2/180*pi;% delta alpha for map searching
d_b=10/2/180*pi;% delta beta for map searching

R2=zeros(nref,1);
cfref=zeros(ns,ns,nref);
%pxzI=zeros(ns,ns,nalpha*nimg,nbeta*ngamma);
pxzI=zeros(ns,ns,nalpha,nref);
pxzC=pxzI;
%pxzIa=zeros(ns,ns,nalpha);
% %%  compute the conjugate for each references in fourier space
% for j=1:nbeta
%     for k=1:ngamma
%         iref=(j-1)*ngamma+k;
%         R=templates(:,:,iref);
%         R2(iref)=R(:)'*R(:); 
%         cfref(:,:,iref)=conj(fftn(templates(:,:,iref)));
%     end;
% end;
%%
Pxz=zeros(ns,ns,nalpha,nbeta*ngamma);

%li=zeros(nalpha*nimg,1);
li=-inf(nalpha*nimg,nref);
lir=zeros(nimg,1);
%lix=lir;liy=lir;
%Firstli=1;
Iclassmean=zeros(ns,ns,nimg,nref);
ICTF=Iclassmean; 
classmean=zeros(ns,ns,nref);
UnrotImgs=zeros(ns,ns,nalpha);
fUnrotImgs=zeros(ns,ns,nalpha);
Px=0;
Pxs=ones(nimg,1);
niter=3; % iteration of estimation
a=1.2*ones(nimg,1);  % ratio estimate, I(i)=a(i)R
SdW=Pxs;
SpxzN=Pxs;
repeat=1; % repeat classmean refinement before reconstruction
Rs=templates;
for iiter=1:niter
  if mod(iiter,repeat)==1 || repeat==1
      % CTF prepare references
    templates=rsMakeTemplates(angleList,vol);
%     for j=1:nbeta
%         for k=1:ngamma
%             iref=(j-1)*ngamma+k;
%             if 1%iiter==1  % is this necessary?
%                 fR=fftn(templates(:,:,iref)).*ifftshift(c(:,:,i));
%                 R=real(ifftn(fR));
%                 templates(:,:,iref)=R;
%                 cfref(:,:,iref)=conj(fR);
%                 R2(iref)=R(:)'*R(:); 
%             else
%                 R=templates(:,:,iref);
%                 R2(iref)=R(:)'*R(:); 
%             end;
%             
%             
%         end;
%     end;
   end;
  A1=zeros(nalpha,nref);
  A2=A1;
  sIR=A1;
  pxzN=A1;
  Px=A1;
  dW3=A1;
    % go through every image

  for i=1:nimg 
          for j=1:nbeta
            for k=1:ngamma
                iref=(j-1)*ngamma+k;
                if 1%iiter==1  % is this necessary?
                    fR=fftn(templates(:,:,iref)).*ifftshift(c(:,:,i));
                    R=real(ifftn(fR));
                    Rs(:,:,iref)=R;
                    cfref(:,:,iref)=conj(fR);
                    R2(iref)=R(:)'*R(:); 
                else
                    R=templates(:,:,iref);
                    R2(iref)=R(:)'*R(:); 
                end;                        
            end;
          end;
          figure(1);
        ImagicDisplay(Rs);
    %if iiter==1  % save fUnrotImgs(i,q) to another stack
        RotImg=rotImgs(:,:,i);  % pick one particle image from the database
        %RrotImg=rsRotateImage(RotImg,);  % rotate one particle image back by the measured alpha
        I2=RotImg(:)'*RotImg(:);    
        %UnrotImgs=rsRotateImage(RrotImg,alphaRs*180/pi); % Rock the rotated image
        UnrotImgs=rsRotateImage(RotImg,ralphas(i)+alphaRs*180/pi); % Rock the rotated image

        for q=1:nalpha        
            fUnrotImgs(:,:,q)=fftn(UnrotImgs(:,:,q));   % compute the fftn for the set of rocking images
        end;
        % pwc=zeros(ns,ns,nalpha,nbeta);
    %end; 
    [pwcf,dWf,Flagf]=P_wc2(alphaRs,betaRs,d_a,d_b,SigmaB,SigmaR,ns,r(i)+dr,y(i),SigmaC);  % it return nalpha*nbeta values
    [pwcr,dWr,Flagr]=P_wc2(alphaRs,betaRs,d_a,d_b,SigmaB,SigmaR,ns,r(i)-dr,-y(i),SigmaC);
    pwc=pwcf*lambda+pwcr*(1-lambda);
    lcz=log(pwc);    
    %lczr=log(pwc*(1-lambda));
    dW2=dWf+dWr;
    
    Flag=(Flagf+Flagr)>0;
    for j=1:nbeta
       
%         imacs(residual>25);
%         title(num2str(j));
%         drawnow;
        if Flag(round(nalpha/2),j)==0  %Flag(1,j)==0 % should test nalpha  when beta=0
            %display([num2str([betaRs(j),betas(i)]*180/pi), ' pwc is too small']);
        else
            
            for k=1:ngamma
                iref=(j-1)*ngamma+k;
%                 if iref==28
%                    display('Here');
%                 end;
                for q=1:nalpha
                    
                    IR=real(fftshift(ifftn(fUnrotImgs(:,:,q).*cfref(:,:,iref))));
                    % does not include pc
                    %lI=-ns*ns*log(sqrt(2*pi)*sigmaN)-(I2+R2((j-1)*ngamma+k)-2*IR)/(2*sigmaN^2);
                    % include pc
                    lI0=(I2+a(i)^2*R2(iref)-2*a(i)*IR);%-ns*ns*log(sqrt(2*pi)*SigmaN)
                    
                    lI=-ns*ns*log(sqrt(2*pi)*SigmaN)-(I2+a(i)^2*R2(iref)-2*a(i)*IR)/(2*SigmaN^2)+lcz(:,:,q,j);
                    maxlI=max(lI(:));
%                     if Firstli  % do this because li is iniated as zeros
%                         li=maxlI+li; % replace the zero matrix with a (probably negative) value
%                         Firstli=0;
%                     else
                        li(nalpha*(i-1)+q,iref)=maxlI;
%                     end;
                    if isinf(maxlI)
                        Pxz(:,:,q,iref)=0;
                    else
                        Pxz1=exp(lI-maxlI); % for each image, each reference
                        Pxz(:,:,q,iref)=Pxz1;
                    end;
                    dW3(q,iref)=dW2(q,j)*sum(Pxz1(:));
                    % now start computing the convolution of image and Pxz
                    % for a set of alphaRs
                    tempf=fUnrotImgs(:,:,q).*conj(fftn(ifftshift(Pxz1))); % weighted and shifted image
                    temp=real(ifftn(tempf.*ifftshift(c(:,:,i)))); % CTF-filtered W&Sed image
                    tempc=ifftshift(c(:,:,i).^2).*conj(fftn(ifftshift(Pxz1)));  % W&Sed CTF
                    tempn=real(ifftshift(lI0.*conj(fftn(ifftshift(Pxz1))))); % can be weighted with easier method
                    pxzI(:,:,q,iref)=temp;
                    pxzC(:,:,q,iref)=fftshift(real(ifftn(tempc)));  % the same as image
                    pxzN(q,iref)=sum(tempn(:));
                    %pxzI(:,:,q,iref)=real(ifftn(tempf.*ifftshift(c(:,:,i))));
                    %pxzC(:,:,q,iref)=tempc(:)'*tempc(:);
                    R=Rs(:,:,iref);                    
                    sIR(q,iref)=R(:)'*temp(:);  % to estimate a(i)
                    %to estimate SigmaC
                    %fftn(R).*conj(fftn(ifftshift(Pxz(:,:,q,iref)))))
                    %pxzIa(:,:,q)
                    
                end;
                
                %pxzIr=sum(pxzIa,3);
            end;
        end;
        
    end;
   
    %[lir(i),lix(i),liy(i)]=max2d(li);  % one lir for each experimental beta/particle/alpha*image lix:alpha, liy:beta, gamma
    lir(i)=max2d(li(nalpha*(i-1)+1:nalpha*(i),:));
    %figure();
    %plot(li);
    
    % compute p(X), the sum of p(X,Z) over Z
     for q=1:nalpha
%         % [lir,lin]=max(li(nalpha*(i-1)+q,:));
%          
          for iref=1:nref
              scaler=exp(li(nalpha*(i-1)+q,iref)-lir(i));
              pxzI(:,:,q,iref)=pxzI(:,:,q,iref)*scaler;
              pxzC(:,:,q,iref)=pxzC(:,:,q,iref)*scaler;
              Pxm=Pxz(:,:,q,iref);
              Px(q,iref)=sum(Pxm(:))*scaler;
              pxzN(q,iref)=pxzN(q,iref)*scaler;
              A1(q,iref)=sIR(q,iref)*scaler; %*Px;              
              A2(q,iref)=R2(iref)*Px(q,iref);  % can be put in the earlier loop
              dW3(q,iref)=dW3(q,iref)*scaler;
              Pxs(i)=Pxs(i)+Px(q,iref);   
          end;
     end;
     a(i)=sum(A1(:))/sum(A2(:));
     
     % normalize the pxzI for each set of images (different alphas?)
%     for iref=1:nref
%         for q=1:nalpha
%             
%             %pxzI(:,:,nalpha*(i-1)+q,iref)=pxzI(:,:,nalpha*(i-1)+q,iref)*exp(li(nalpha*(i-1)+q,iref)-liM)/(exp(lir(i)-liM)*Pxs(i));
%             pxzI(:,:,q,iref)=pxzI(:,:,q,iref)*exp(li(nalpha*(i-1)+q,iref)-lir(i))/Pxs(i);
%             %dW2(q,iref)=dW2(q,iref)/(2Pxs(i)*2*N);
%         end;
%     end;
    SdW(i)=sum(dW3(:));
    SpxzN(i)=sum(pxzN(:));
    pxzI=pxzI/Pxs(i);
    pxzC=pxzC/Pxs(i);
    Iclassmean(:,:,i,:)=squeeze(sum(real(pxzI),3));
    ICTF(:,:,i,:)=squeeze(sum(real(pxzC),3));
  end;
  lim=max(lir); 
  Pxsm=Pxs.*exp(lir-lim);
  SdW=SdW.*exp(lir-lim);
  SigmaC=sqrt(sum(SdW(:))/(2*sum(Pxsm(:))*nimg))%
  SigmaN=sqrt(sum(SpxzN.*exp((lir-lim))/sum(Pxsm(:)))/(nimg*ns*ns))
  SigmaN=1;
%[liM,lx,ly]=max2d(li)
%a
%% get the classmean
%classmean=squeeze(sum(real(pxzI),3));
classmean=squeeze(sum(Iclassmean,3));
norms=squeeze(sum(ICTF,3));
DrawClassmean(classmean,angleList,li-lim);
mbw=ImageArray(classmean);
figure(3);
subplot(111);
imacs(mbw);
SetGrayscale;
drawnow;
if mod(iiter,repeat)
    display('Refine classmeans');
else
    display('Make a new model and projections');   
    [vol norm]=rsGriddingReconstruction(angleList,classmean,norms,k);
    % vol=rsFourierReconstruction(angleList, classmean, CTF,symmetry,kWiener);
    figure(4);
    ShowSections(vol);
    figure(5);
    fnorm=fftshift(abs(fftn(norm)));
    ShowSections(fnorm);
    drawnow;
    %templates=rsMakeTemplates(anglesList,vol);    
end;
%DrawClassmean(classmean,angleList,lir); 
%DrawPxz(Pxz,li-max(lir));
    %Q=log(Pxz)*P

% class mean will be convolution of Pzx and image

end;