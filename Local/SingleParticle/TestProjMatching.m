% TestProjMatching

load ~/Desktop/ministack
% st0=rawStack(:,:,1:2000);
st0=rawStack;

[sp st1]=spNormalizeAndCenter(st0,5.8);
st1=spDownsample(st1,24);
%%
vol=Downsample(fuzzymask(32,3,[11 6 4],4),24);
save vol vol

volmsk=fuzzymask(24,3,[10 7 6],2);
% ShowSections(vol);
angStep=15;

 iter=1;

%% Refinement loop



while iter<10
    load vol
    iter
    
    [projs angles]=spReproject(vol.*volmsk,angStep);
    
    numProjections=size(angles,1)
    figure(1);
    imats(projs,2);
    drawnow;
    
    %% Rotational alignment of images
    tic
    [sp aliStack]=spAlignRot(st1,sp,projs,1);
    toc
    %%
    [classes ncls]=spSumClasses(aliStack,sp);
    figure(2);
    imats(classes,2);
    drawnow;
    
    vol=spReconstruct(classes,ncls,angles,.001);
%      vol=GaussFilt(vol,.2);
    figure(3);
    ShowSections(vol);
    save vol vol;
    
    iter=iter+1;
end;