% rsShowTemplateExamples.m
%  Give examples of projections of a 3D model at given beta and gamma
%  angles, no CTF, displayed in a resizable window.

pixA=5.8;  % desired pixel size
resolution=30;  % angstroms
outName='~/Desktop/templates';

mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';

    
    % Load the 3D map
    disp('Loading the 3D map');
    [origMap s]=ReadMRC(mapName);
    mpixA=s.pixA;
    nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
    nt=ceil(nt1/8)*8;
    [map finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
    map=map*pixA;  % approx amplitude correction (V-A scaling)

    betas=[45 90 105 120 135];
    gammas=[0 30 60 90 120 150];  % 2-fold symmetry : gamma<180
    nb=numel(betas);
    ng=numel(gammas);
    angles=zeros(nb*ng,3);
    for i=1:nb
        angles((i-1)*ng+1:i*ng,2)=betas(i);
        angles((i-1)*ng+1:i*ng,3)=gammas;
    end;
    angles

    % allTemplates=rsMakeTemplatesQuick(angleList,map);
    templates=rsMakeTemplates(angles,map);
    fc=pixA/resolution;
    filtTemplates=SharpFilt(templates,fc,fc/10);
    
    mag=1;  % scale of display.  1 means 1:1 pixel
    ImagicDisplay1(filtTemplates,mag);
    
    WriteImagic(filtTemplates,outName);