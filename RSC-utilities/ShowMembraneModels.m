% ShowMembraneModels.m

% maxNFiles=200;
%
%     rootPath=uigetdir('.','Select a directory containing mi files');
%     if isnumeric(rootPath)
%         return
%     end;
%     cd(rootPath);
%
% %     names=FindFilenames('./','.+mi\.mat');
%     names=FindFilenames('./','.+mi\.txt');
maxNFiles=2000;

[names,pa]=uigetfile('*mi.txt','multiselect','on');
[rootPath,infoPath]=ParsePath(pa);
cd(rootPath);
%%
models=[];
amps=[];
vesRs=[];
vesAs=[];
vesOks=false(1,4);
j=1;
%     char(names)
for i=1:min(maxNFiles, numel(names));
    mi=ReadMiFile([infoPath names{i}]);
    plot(mi.vesicleModel);
    title([num2str(i) '  ' names{i}],'interpreter','none');
    drawnow;
    
    models(:,i)=mi.vesicleModel;
    amps(i)=mi.ctf(2).ampFactor;
    nves=numel(mi.vesicle.x);
    if nves>0
        vesRs(j:j+nves-1,1)=real(mi.vesicle.r(:,1));
        vesAs(j:j+nves-1,1)=real(mi.vesicle.s(:,1));
        vesOks(j:j+nves-1,:)=mi.vesicle.ok;
        j=j+nves;
    end;
    %         drawnow;
end;
%%
plot(models);
meanModel=mean(models,2);
plot(meanModel);
return

%%     Write models into every info file
minAmp=.0035;
maxAmp=.006;
minR=100/mi.pixA
maxR=300/mi.pixA

d=dir(infoPath);
for i=1:numel(d)
    name=d(i).name;
    if strndcmp(name,'mi.txt')
        disp(name);
        mi=ReadMiFile([infoPath name]);
        nves=numel(mi.vesicle.x);
        if nves>0
            amps=real(mi.vesicle.s(:,1));
            rads=real(mi.vesicle.r(:,1));
            flags=amps>minAmp & amps <maxAmp & ...
                rads>minR & rads<maxR & mi.vesicle.ok(:,1);
            mi.vesicle.ok(:,2)=flags;
            mi.vesicle.ok(:,3:4)=true;
            
            %             mi.vesicleModel=meanModel;
            %             mi.ctf(2).ampFactor=1.14;
round([amps*1e4 rads mi.vesicle.ok])
           WriteMiFile(mi,[infoPath name]);
        end;
    end;
end;
