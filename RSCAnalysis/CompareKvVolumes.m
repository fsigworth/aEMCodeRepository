% reCompareKvVolumes

 cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j4');
 [v5,s5]=ReadMRC('mrc/i40av01.mrc');
 
 cd('/Users/fred/EMWork/Hideki/140625n/Reconstructions/Recon112j6')
 [v4,s4]=ReadMRC('mrc/i40av01.mrc');
 
 load('/Users/fred/aEMCodeRepository/AMPAR/KvMap.mat');
 s0.pixA=pixA;
%%
    theta=-5;
    
v0=rsRotateImage(Crop(map,128),50);
sigma0=.5;
mag=.53;
 sz=size(v4);
 
 v1=GaussFilt(DownsampleGeneral(v0,72,mag),.2);
 v1s=circshift(v1,[0 0 0])+sigma0*randn(sz);
 ct=37;
 slice=[ct ct 26];
 
 figure(11);
 ShowSections(v1s,slice,45);
 drawnow;
%  
 v5s=MirrorX(rsRotateImage(v5,theta));
 v5s=5*DownsampleGeneral(v5s-.8*GaussFilt(v5s,.2),72,1.0);
 figure(12);
 ShowSections(v5s,slice,45);
 drawnow;
 
 figure(13);
 ShowSections(v5s-.28*v1s,slice,45);
 title(theta);
 drawnow;
 
 %%
 
 v4s=v4+GaussHP(v4,.2);
 v1u=Downsample(v1s,112);
 v4u=Downsample(v4s,112);
 v5u=Downsample(v5s,112);
 
 zlim=28:100;
 
 figure(14);
 ImagicDisplay3(v1u(:,:,zlim),0,1);
 
 figure(15);
 ImagicDisplay3(v4u(:,:,zlim),0,1);
 
 figure(16);
  ImagicDisplay3(v5u(:,:,zlim-1),0,1);

 %%
%  Write these into the 140625n folder.
WriteMRC(v1s,s4.pixA,'v0.mrc' );
WriteMRC(v4s,s4.pixA,'v4.mrc');
WriteMRC(v5s,s5.pixA,'v5.mrc');
