% SmileyCrossCorr.m

vlSetDarkGraphics(18,[.28 .28 .3]);
vlSet1080Figure

load smiley256.mat % variable smx

figure(1);

n=512;
ns=64;
sm=1-DownsampleGeneral(GaussFilt(smx,.1),ns,1/2.5);
imags(sm);

ref=Crop(sm,n);
pos=[100 400; 150 150];
img=0;
for i=1:size(pos,1);
    img=img+ExtractImage(sm,pos(i,:),n,1);
end;

mysubplot(131);
imags(ref);
axis off;

mysubplot(132);
imags(img)
axis off;

cc=real(ifftn(fftn(img).*conj(fftn(ifftshift(ref)))));
mysubplot(133);
imags(cc);
axis off;
return

%%
imgn=img+2.5*randn(n);
mysubplot(132);
imags(GaussFilt(imgn,.2))
axis off;

ccn=real(ifftn(fftn(imgn).*conj(fftn(ifftshift(ref)))));

mysubplot(133);
imags(ccn);
axis off;

return

%% Einstein from noise
makeVideo=1;
vName='Einstein';
if makeVideo
    v=VideoWriter([vName '.mp4'],'MPEG-4');
    v.FrameRate=15;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;



fref=conj(fftn(ifftshift(ref)));
acc=zeros(ns);
acc=zeros(n);
ctr=ceil((n+1)/2);
step=1;
for i=1:2000
    if i==201
        step=10;
    end;
    img=randn(n);
    cc=real(ifftn(fftn(img).*fref));
    [val,ix,jy]=max2d(cc);
    imgsh=circshift(img,[ctr-ix ctr-jy]);
    acc=acc+imgsh;
    %     acc=acc+ExtractImage(img,[i j],ns);
    if mod(i,step)==0
        mysubplot(133);
        imags(acc);
        axis off;
        text(round(.7*n),round(.9*n),[num2str(i) ' summed'],'color','w','fontsize',16);
        mysubplot(132);
        imags(img);
        hold on;
        plot(ix,jy,'y+');
        hold off;
        axis off;
        drawnow;
        if makeVideo
            f=getframe(gcf);
            writeVideo(v,f);
        end;
        
    end;
end;
    close(v);


