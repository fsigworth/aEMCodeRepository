% rlLoadStarImages

findClasses=[];

disp('Get a star file');
[nm,pa]=uigetfile('*.star');
disp(['Reading ' pa nm]);
[names,data]=ReadStarFile([pa,nm]);
% rootPath=ParsePath(ParsePath(pa));
% cd(rootPath)

if numel(strfind(names{1},'opt'))>0 % we have an optics block
    d=data{2}; % 2nd block will contain particle or micrograph data.
else
    d=data{1};
end;
nl=numel)(d.rlnMicrographName);

[uNames,uLines]=unique(d.rlnMicrographName);
nu=numel(uNames);

for i=1:nu
    if exist(uNames{i},'file')
        mic=ReadMRC(uNames{i});
        
    cPars=rlCTFParsFromStar3Line(names,data,uLines(i),bValue);
    


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

