% MakeDistribution.m

disp(' ');
disp('MakeDistribution.m');

sourcePath='~/aEMCodeRepository/';
targetPath='~/aEMCodeDistribution/';

CheckAndMakeDir(targetPath);

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
         };
     
     system(['cp ' sourcePath 'startupDistribution.m ' targetPath 'startup.m']);
     system(['mkdir ' targetPath 'aLibs']);
     
     for i=1:numel(folders)
         name=folders{i};
         str=['cp -R ' sourcePath name ' ' targetPath name];
         disp(str);
         system(str);
     end;
cd ~/
str='tar -czf aEMCodeDistribution.tgz aEMCodeDistribution';
disp(str);
system(str);
disp('done.');