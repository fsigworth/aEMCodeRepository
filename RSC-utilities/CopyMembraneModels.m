% CopyMembraneModels.m

%BValue=40;

disp('Get the reference membrane model');
[refName,pa]=uigetfile('*mi.txt');
if isnumeric(refName)
    return
end;
[rootPath,infoPath]=ParsePath(pa);
cd(rootPath);
%
miRef=ReadMiFile([infoPath refName]);

frameSets=miRef.frameSets

model=miRef.vesicleModel;
figure(1); clf;
plot(model);

disp('Get the first file to update');
[firstName,pa]=uigetfile('*mi.txt');
if isnumeric(firstName)
    return
end;
[rootPath,infoPath]=ParsePath(pa);
cd(rootPath);


% plot(model);
% title([refName],'interpreter','none');
% drawnow;
startInd=1;
%%
    d=dir(infoPath);
    for i=startInd:numel(d)
        name=d(i).name;
        if strndcmp(name,'mi.txt')
            mi=ReadMiFile([infoPath name]);
            if size(mi.vesicle.s,2)<2 % not refined
                %             if ~isfield(mi,'frameSets') || mi.frameSets(1,2)<15
                %                 mi.frameSets=frameSets; disp([name ' changed']);
                mi.vesicleModel=model;
                disp(name);
                %             nCtf=numel(mi.ctf);
                % %             for i=1:nCTF
                %                 mi.ctf(i).B=BValue; mi.ampFactor=1;
                WriteMiFile(mi,[infoPath name]);
                disp(name);
            end;
        end;
    end;

