% startup.m for Farnam
disp('Matlab startup');
% basePath='/ysm-gpfs/datasets/cryoem/cryoMatlab/';
basePath='~/Matlab_EM/';
folders={'aLibs/EMBase'
         'aLibs/EMCCD'
         'aLibs/EMIO'
         'aLibs/EMIOUtils'
         'aLibs/EMSpec'
         'aLibs/GriddingLib'
         'aLibs/MEX'
         'AMPAR'
         'FourierReconstruction'
         'K2Camera'
         'Krios'
         'MultiExposure'
         'Relion'
         'Realtime'
         'RSC'
         'RSC-Apps'
         'RSC-utilities'
         'RSCAnalysis'
         'RSCPicking'
         'RSCReconstruct'
         'RSCSimulate'
         'RSCVesicleFinder'
         };
         
for i=1:numel(folders)
    addpath([basePath folders{i}]);
end;
