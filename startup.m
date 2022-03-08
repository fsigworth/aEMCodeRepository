% startup.m
% Put this script somewhere where Matlab will look on startup, e.g. in home
% folder or in ~/Documents/MATLAB
host='MyHost';
homePath='~/';
basePath='~/aEMCodeRepository/'; % We asssume that the code is here.

disp([which('startup') ' on ' host]);

setenv('HOSTNAME',host);
cd(homePath);
folders={['aLibs' filesep 'EMBase']
         ['aLibs' filesep 'EMCCD']
         ['aLibs' filesep 'EMIO']
         ['aLibs' filesep 'EMIOUtils']
         ['aLibs' filesep 'EMSpec']
         ['aLibs' filesep 'GriddingLib']
         ['aLibs' filesep 'MEX'   ]   
         ['aLibs' filesep 'Others']
         'EMClass'
         'FourierReconstruction'
         'K2Camera'
         'MultiExposure'
         'Relion'
         'Realtime'
         'RSC'
         'RSC-Apps'
         'RSC-utilities'
         'RSCPicking'
         'RSCReconstruct'
         'RSCSimulate'
         'RSCVesicleFinder'
         'RSCAnalysis'
         'AMPAR'
         'Kv'
         'VesicleML'
%         ['aLibs' filesep 'Others' filesep 'bfmatlab']
%          'RSCAdaptiveBasis'
%          'RSCAdaptiveBasis/utils'
         };
         
for i=1:numel(folders)
    addpath([basePath folders{i}]);
end;

