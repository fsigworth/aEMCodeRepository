% FT2SmileyRot.m
% Example of rotation propert of FT
%
vlSetDarkGraphics(20);
clf;
vlSet1080Figure(1,1,[1600 640]);

rotStep=4; % degrees
nRots=round(360/rotStep);
n=128;
makeVideo=1;


% cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.4FT2D/Figs/')
cd('/Users/fred/Documents/Documents - Katz/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.4FT2D/Figs/')
vName='FT2SmileyRot.mp4';
if makeVideo
    v=VideoWriter(vName,'MPEG-4');
    v.FrameRate=10;
    open(v);
    disp(['Making movie file: ' vName '.mpr']);
end;

load smiley256.mat % loads smx

sm=Downsample(Crop(-smx,200),n);
ctr=ceil((n+1)/2);

smf=fftshift(fftn(ifftshift(sm)));
smf(ctr,ctr)=0;

subplot(121);
imags(sm);
axis equal;
subplot(122);
imacx2(smf,.3);
axis equal
msk=fuzzymask(n,2,.48*n,.1);


for i=0:nRots
    smr=msk.*rsRotateImage(sm+1,i*rotStep);
    smf=msk.*fftshift(fftn(ifftshift(smr)));
    smf(ctr,ctr)=0;

    subplot(121);
    imags(smr);
    axis equal;
    subplot(122);
    imacx2(smf,.3);
    axis equal
    drawnow;

    if makeVideo
        f=getframe(gcf);
        writeVideo(v,f);
    else
        pause;
    end;
end;

if makeVideo
    close(v);
    disp(['Movie ' vName ' written.'])
    disp(['In: ' pwd])
end;



