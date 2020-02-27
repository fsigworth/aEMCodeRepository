% startup.m
host='Katz';
homePath='~/';
basePath='~/aEMCodeRepository/';

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
         ['aLibs' filesep 'Others' filesep 'bfmatlab']
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
%          'RSCAdaptiveBasis'
%          'RSCAdaptiveBasis/utils'
         };
         
for i=1:numel(folders)
    addpath([basePath folders{i}]);
end;
% disp('Using desktop RSCReconstruct');
% addpath('/Users/fred/Desktop/RSCReconstruct');
% addpath('/Users/fred/Dropbox/OHSU_images/matlabprogram');
% disp('Using Liguo''s code');

