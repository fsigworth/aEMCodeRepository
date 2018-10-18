% SortComparedFiles.m
% From the inds.mat indices written by CompareJpegs, sort the mi and merged
% images into a common directory.

basePath='/ysm-gpfs/scratch60/fjs2/160909/';
findVesicleDirs={'KvLipo121_2w10_v3.5/';'KvLipo121_2w10_v3.6/';
    'KvLipo121_2w10_v3.7/';'KvLipo121_2w10_v3.8/'};
outputDir='KvLipo121_2w10v3t/';
cd(basePath);
CheckAndMakeDir(outputDir,1);
CheckAndMakeDir([outputDir 'Info/']);
CheckAndMakeDir([outputDir 'Merged/']);

load inds.mat % get inds and names
nf=numel(names);
ni=numel(inds);
if ni<nf
    inds(nf)=0;
end;

% set defaults
inds=max(1,inds);
ndir=numel(findVesicleDirs);

for i=1:nf
    j=inds(i);
    if j<=ndir
        jname=names{i};
        baseName=jname(1:end-6);
        miName=[baseName 'mi*'];
%         if exist(miName,'file')
%             mi=ReadMiFile(miName);
            str1=['cp -v ' findVesicleDirs{j} 'Info/' miName ' ' outputDir 'Info/'];
            str2=['cp -v ' findVesicleDirs{j} 'Merged/' baseName 'm* ' basePath outputDir 'Merged/'];
%              disp(str1);
%             disp(str2);
            system(str1);
            system(str2);
%         else
%             disp(['Not found: ' miName]);
%         end;
        
    end;
end;




% cd([basePath findVesicleDirs{1}]);
% d=dir('Jpeg');
% names=cell(0,1);
% for i=3:numel(d)
%     names{i-2}=d(i).name;
% end;
% nf=numel(names);
