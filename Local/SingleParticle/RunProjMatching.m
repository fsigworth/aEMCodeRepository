% RunProjMatching

mode=2;

% load ~/Desktop/ministack  % 8000 particles
load ~/Desktop/stackDSR32  % 22000 particles

[sp st1]=spNormalizeAndCenter(rawStack,5.8);
st1=GaussHP(st1,.05,1);
figure(1);
imats(st1,2);
drawnow;

switch mode
    case 1
        %%  Initial volume model is an ellipsoid
        
        vol=fuzzymask(32,3,[11 6 4],4);
        save vol vol
        
        % mask for the volume
        volmsk=fuzzymask(32,3,[14 10 8],2);
        
    case 2
        %%  Initial volume model is from pdb
        load ~/Desktop/DSRmap32  % 32^3 map
        
        vol=SharpFilt(sim2,5.8/50,.02);  % filter to 40 Å
        ShowSections(vol);
        % binary mask
        volmsk=vol>120;
        
        WriteMRC(vol,5.8,'MyVolumeOrig.mrc');
end;

save vol vol


angStep=10;  % angular step, in degrees

iter=1;

%% Refinement loop


while iter<20
%     load vol
    iter
    
    [projs angles]=spReproject(vol.*volmsk,angStep);
    
    numProjections=size(angles,1)
    figure(1);
    imats(projs);
    drawnow;
    
    %% Rotational alignment of images.
    %     This is primitive because to save time we don't do translational
    %     alignment.
%     [sp aliStack]=spAlignRot(st1,sp,projs,1);
    [sp aliStack]=spAlignRotTrans(st1,sp,projs);
%%
% [aliStack dxy]=spMatchTrans(aliStack,sp,projs);  % match displacements
% %%
% figure(1);
%     clf;
%     plot(dxy(:,1),dxy(:,2));
%     drawnow;

%%
%  [sp aliStack]=spAlignRot(aliStack,sp,projs,1);

    %%
    [classes ncls]=spSumClasses(aliStack,sp);
    figure(2);
    imats(classes);
    drawnow;
    
    vol=spReconstruct(classes,ncls,angles,.001);
        vol=GaussFilt(vol,.25);
    figure(3);
    ShowSections(vol);
    save vol vol;  % save this in case we want to restart.
    
    WriteMRC(vol,sp.pixA,['MyVolumeRT' num2str(iter) '.mrc']);  % save for viewing with Chimera
    
    iter=iter+1;
end;

