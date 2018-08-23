% We use the origImg and mxCC8 cross-correlation image, along with the mi
% structure.  mxCC8 is downsampled by 2 from the origImg.

cc=single(mxCC8)/100;
dimg=imscale(GaussHP(GaussFilt(origImg,.2),.005),256,.01);
% dimg=imscale(GaussHP(GaussFilt(mVesicles+mParticles,.2),.005),256,.01);
dimgs=imscale(GaussHP(GaussFilt(origImg-meMakeModelVesicles(mi,2048),.2),.005),256,.001);
dcc=imscale(cc,256,.001);

figure(1); clf; SetGrayscale;
imac(dimg);
%%  Crop out an interesting region

x0=850;
y0=100;
% y0=0;
dx=512;

dimgx=dimg(x0+1:x0+dx,y0+1:y0+dx);
dimgsx=dimgs(x0+1:x0+dx,y0+1:y0+dx);
dccx=dcc(x0/2+1:x0/2+dx/2,y0/2+1:y0/2+dx/2);
ccx=cc(x0/2+1:x0/2+dx/2,y0/2+1:y0/2+dx/2);
%  cgx=mxGamma(x0/2+1:x0/2+dx/2,y0/2+1:y0/2+dx/2);
cax=mxTempInds(x0/2+1:x0/2+dx/2,y0/2+1:y0/2+dx/2);

% figure(3);
% imac(dccx);
figure(1);
imac(dimgx);

% %% Make a contour plot of the region
% figure(2);
% colormap jet;
% lccx=ccx;
% lccx(lccx<.2)=.2;
% lccx=GaussFilt(lccx,.2);
% % lccx(ccx>1)=1;
% contourf(rot90(lccx',0));
% colorbar
dimgy=dimgx;
dimgsy=dimgsx;
ccy=ccx;
ccf=ccy;
dx1=dx/2;
cay=cax;

%%  Now crop to an even smaller region
x0=0;
y0=128;
% y0=0;
dx=384;
x1=x0/2;
y1=y0/2;
dx1=dx/2;

ccy=ccx(x1+1:x1+dx1,y1+1:y1+dx1);
figure(4);
SetGrayscale;
imacs(ccy);

% cgy=cgx(x1+1:x1+dx1,y1+1:y1+dx1);
cay=cax(x1+1:x1+dx1,y1+1:y1+dx1);
% figure(5);
% imac(cgy);
% colorbar;

% Make a contour plot of the CCF
dimgy=dimgx(x0+1:x0+dx,y0+1:y0+dx);
dimgsy=dimgsx(x0+1:x0+dx,y0+1:y0+dx);
ccf=GaussFilt(ccy,.2);
ccf=ccy;
figure(2);
lccf=ccf;
lccf=GaussFilt(ccf,.2);
lccf(ccf<.4)=.2;
contourf(rot90(lccf',0));
% imacs(lccf);
colorbar
%%  Do the picking on this small region
% assumes these images:
% dimgy, dimgsy, ccf

compImg=0*dimgy;
ccf1=ccf;
di=zeros(dx,dx,3);
figure(6); SetGrayscale;
s0=20;
square=zeros(s0,s0);
square(1,:)=1; square(s0,:)=1;
square(:,1)=1; square(:,s0)=1;
msk=compImg;
mskr=7;
scl=4;
for k=1:50
    
    
    [v i j]=max2d(ccf1);
    ccf1=ccf1.*(1-fuzzymask(dx1,2,mskr,1,[i j]));
%     figure out which reference matched the best
%     g=cgy(i,j);
%     gamma=mod(g,nGamma);
%     h=floor((g-1)/nGamma);
%     hemi=mod(h,2)+1;  % 1 means upper hemisphere, 2 means lower
    ang=cay(i,j);
    ang
    bestTemplate=xTemplates(:,:,ang);
    
    
    
    ptcl=Downsample(xTemplates(:,:,ang),nt*2);
    point=[i j]*2-1;
%     if (v<0.9 && v>.6)  % accepted range of CC peaks
    if (v<1.1 && v>.8)  % accepted range of CC peaks
%         g
        compImg=compImg-ExtractImage(scl*ptcl,2*[i j+10]-1,dx,1);
        compImg=0*compImg;  % Suppress drawing the particle above the box
        msk=msk | (compImg<-4);
        msk=Mask(msk,point,1+0*square,square);
        subplot(2,2,1);
        di(:,:,1)=rot90(dimgy,1);
        di(:,:,2)=rot90(dimgy.*(1-msk)+compImg.*msk,1);
        di(:,:,3)=rot90(dimgy.*(1-msk)+compImg.*msk,1);
        image(uint8(di));
        
        subplot(2,2,2);
        di(:,:,1)=rot90(dimgsy,1);
        di(:,:,2)=rot90(dimgsy.*(1-msk)+compImg.*msk,1);
        di(:,:,3)=rot90(dimgsy.*(1-msk)+compImg.*msk,1);
        image(uint8(di));
        % imacs(compImg); %%
        title(v)
        subplot(2,2,3);
imacs(ccf1);
        
        drawnow;
    end;
end;
subplot(2,2,3);
imacs(ccf);
subplot(2,2,4);
imacs((ccf>.8) + (ccf>1.2))


minv=v