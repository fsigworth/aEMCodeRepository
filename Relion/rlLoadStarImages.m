% rlLoadStarImages

findClasses=[];

% disp('Get a star file');
% [nm,pa]=uigetfile('*.star');
% disp(['Reading ' pa nm]);
% [names,data]=ReadStarFile([pa,nm]);
% rootPath=ParsePath(ParsePath(pa));
% cd(rootPath);
s=data{1};
fieldNames=fieldnames(s);
nl=numel(s.(fieldNames{1}));
disp([num2str(nl) ' lines']);

disp('Scanning names');
inds=zeros(nl,1);
names=cell(nl,1);
for i=1:nl
    [inds(i),names{i}]=rlDecodeImageName(s.rlnImageName{i});
end;
[uniqueFileNames,[],fi2]=unique(names);
nuf=numel(uniqueFileNames);
plot(inds)

% We're in a Stack directory. Look for image files in equivalent
% directories
ourDir=pwd;
fileFound=false(nuf,1);
cd ../..
d=dir;
for i=3:numel(d)
    if d(i).isdir


names{10000}

