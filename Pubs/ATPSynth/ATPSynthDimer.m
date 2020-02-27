% ATPSynthDimer.m


cd('/Users/fred/EMWork/Nelli/ATPSynthDimer20180315/') % base directory

% Our experimental map
[me,se]=ReadMRC('20180315/Stack2/job022/postprocess.mrc'); % experimental map
ne=size(me,1); % 128
nex=1.5*ne; %192
nesm=nex/2; % small cropped size = 96
mex=Crop(me,nex);
mesm=Crop(me,nesm);

figure(1);
ShowSections(me);
drawnow;


%% Dimer map
[md,sd]=ReadMRC('emd_8151.map');
nd=size(md,1);
mds=DownsampleGeneral(GaussFilt(md,.1*sd.pixA),ne,sd.pixA/se.pixA);
figure(2);
ShowSections(mds);
drawnow;


%% Make a small map from the pdb file.

% mes=DownsampleGeneral(me,nd,se.pixA/sd.pixA);
fc=.05; % Gaussian lowpass for atomic model
[composite,protein]=SolventAndProteinDensity('DimerFigure/5ARA.pdb');
refsm=GaussFilt(DownsampleGeneral(composite-composite(1),ne,1/se.pixA),fc*se.pixA);

figure(3);
ShowSections(refsm);

%%
% Write out 128^3 maps, pixA=4.1
WriteMRC(me,se.pixA,'ExpMonoMap.mrc');
WriteMRC(refsm,se.pixA,'PdbMonoMap.mrc');
WriteMRC(mds,se.pixA,'DimerMap.mrc');
