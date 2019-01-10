function k2PipelineDisplay

pa.ignoreCTF=0;
archiveDir='tempGraphics/';
deleteOriginals=0;

% We assume we're in a proper directory
disp('k2PipelineDisplay working directory is');
disp(pwd);
d=dir;
figure(1);
firstRun=1;
while true
    for i=1:numel(d)
        dirName=d(i).name;
        switch lower(dirName)
            case 'jpeg'
                arch=[AddSlash(dirName) archiveDir];
                if ~deleteOriginals
                    CheckAndMakeDir(arch,1);
                end;
                ScanForDisDats(dirName,firstRun,pa,deleteOriginals,...
                    arch);
                firstRun=0;
        end;
    end;
    pause(5);
end
end

function ScanForDisDats(dirName,firstRun,pa,del,archDir)
d=dir(dirName);
for i=1:numel(d)
    name=d(i).name;
    p=strfind(name,'DisDat.mat');
    if numel(p)>0 && ~(pa.ignoreCTF && strndcmp(name,'ctfDisDat.mat',13))
        datName=[AddSlash(dirName) name];
        try
            s=load(datName);
            
            names=fieldnames(s);
            for j=1:numel(names)
                dName=names{j};
                clf;
                DrawFigureFromData(s.(dName));
                if j==1
                    jpegName=name(1:p(end)-1);
                else
                    jpegName=[name(1:p(end)-1) num2str(j)];
                end;
                jpegNamex=[AddSlash(dirName) jpegName '.jpg'];
                print('-djpeg','-r150',jpegNamex);  % save the CTF window.
            end;
            disp(['Converted to jpg: ' datName]);
            if del
                delete(datName);
                disp('...deleted');
            else
                str=['mv ' datName ' ' archDir name];
                disp(str);
                system(str);
            end;
            firstRun=0;
        catch
            disp(['waiting for file: ' datName]);
        end;  % try
    elseif firstRun
       % disp('No DisDat.mat files found on first scan.');
    end;
end;
end

