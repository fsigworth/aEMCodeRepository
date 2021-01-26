% rlMicrographTrends
% 
% Try to determine bad micrographs


CTFDir='CTFFind/job026/';
[nms,dat]=ReadStarFile([CTFDir 'micrographs_ctf.star']);

% Get allMis
load Info_7_picked_all/allMis.mat

%%
opt=dat{1};
mic=dat{2};

figure(1)
subplot(311);
plot(mic.rlnCtfFigureOfMerit);
subplot(312);
plot(mic.rlnCtfMaxResolution);
subplot(313);
plot(mic.rlnDefocusU/1e4);

figure(2);
hist(mic.rlnCtfMaxResolution,1000);

mc=struct;
mc.names=mic.rlnMicrographName;
mc.defs=(mic.rlnDefocusU+mic.rlnDefocusV)/2e4;
mc.res=mic.rlnCtfMaxResolution;
mc.fom=mic.rlnCtfFigureOfMerit;

nl=numel(mc.names);
mc.baseNames=cell(nl,1);
nErrs=0;
for i=1:nl
    [~,mc.baseNames{i}]=fileparts(mc.names{i});
%     If there are no errors, the mi and star names match!
    if ~strcmp(mc.baseNames{i},allMis{i}.baseFilename)
        nErrs=nErrs+1;
        if nErrs<10
            disp(['Error at ' numtstr(i)]);
        end;
    end;
end;

% Insert extra CTF fields
nmi=numel(allMis);
for i=1:nmi
    ctf=allMis{i}.ctf(1);
    ctf.resLimit=mc.res(i);
    ctf.ccc=mc.fom(i);
    allMis{i}.ctf=ctf;
end;
    
%%

