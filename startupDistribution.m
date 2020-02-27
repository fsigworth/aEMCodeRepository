% startup.m

% Variables for user to set:
defaultHostname='Katz';
homePath='~/';
basePath='~/aEMCodeDistribution/';


host=getenv('HOSTNAME');
if numel(host)<2
    host=defaultHostname;
    setenv('HOSTNAME',host);
end;
disp([which('startup') ' on ' host]);

cd(homePath);
folders={['aLibs' filesep 'EMBase']
         ['aLibs' filesep 'EMCCD']
         ['aLibs' filesep 'EMIO']
         ['aLibs' filesep 'EMIOUtils']
         ['aLibs' filesep 'EMSpec']
         ['aLibs' filesep 'GriddingLib']
         ['aLibs' filesep 'MEX'   ]   
%          ['aLibs' filesep 'Others']
%          ['aLibs' filesep 'Others' filesep 'bfmatlab']
%          'FourierReconstruction'
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
%          'AMPAR'
%          'RSCAdaptiveBasis'
%          'RSCAdaptiveBasis/utils'
         };
         
for i=1:numel(folders)
    addpath([basePath folders{i}]);
end;
% disp('Using desktop RSCReconstruct');
% addpath('/Users/fred/Desktop/RSCReconstruct');
% addpath('/Users/fred/Dropbox/OHSU_images/matlabprogram');
% disp('Using Liguo''s code')
disp('done.');
return

