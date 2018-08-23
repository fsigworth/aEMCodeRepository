% ScanMiParameter.m

% field='frameDose';
field='Defocus, µm';
lastFieldName='noiseModelPars';
% field='vesicle.r';
nv=2;  % pick up this many values

% Have the user select some mi files
[fname, pa]=uigetfile({'*mi.txt' '*mi.mat'},'Select mi files to examine','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
%%
cd(pa);
nim=numel(fname);
val=zeros(nim,nv);
val2=zeros(nim,10);
pt=zeros(nim,1);
mis=cell(nim,1);
for i=1:nim
    %     load(fname{i});
    mi=ReadMiFile(fname{i},lastFieldName);
    mis{i}=mi;
    if isfield(mi,'frameSets') && size(mi.frameSets,2)>1
        pt(i)=mi.frameSets(1,2);
    end;
        if isfield(mi,'ctf') && isfield(mi.ctf,'defocus')
        for j=1:numel(mi.ctf)
            val(i,j)=mi.ctf(j).defocus;
        end;
        end;
    disp([fname{i} '  ' num2str(pt(i),2) '  ' num2str([val(i,:) diff(val(i,:))],'%8.2f')]);
end;
save AllMis.mat mis

%%
valx=val;
if size(val,2)>1
    valx=[val diff(val,1,2)];
end;
subplot(2,1,1);
plot(valx);
ylabel(field);
xlabel('Micrograph');
subplot(2,1,2);
hist(val,70);
xlabel(field);
ylabel('Frequency');

return


%% handling of the 161101/KvLipo122_4b infos:
q=[val diff(val,1,2)];
figure(5)
plot(q);
figure(6)
hist(q(:,3),100);

% Check the defocus jump size
ok=q(:,3)>6.3;
ok=ok & q(:,3)<6.5;

figure(7);
plot(pt);
okp=(pt>39) & (pt<43);

qx=q;
qx(~(okp & ok),:)=NaN;
figure(5);
plot(qx);
okd=(qx(:,1)>1.7) & (qx(:,1)<4);
oka=okp & ok & okd;
sum(oka)
numel(oka)
%%
doMove=0;
for i=1:nim
    name=fname{i};
    if ~oka(i)  % bad file, let's move it.
        [pa,nm,ex]=fileparts(name);
        if strcmp(ex,'.txt')  % currently marked good
            iName=[nm '.txti'];
            str=['mv ' name ' ' iName];
            disp(str);
            if doMove
                system(str);
            end;
            
        end;
    end;
end;


% %%
% b=0;
% while b~='q'
%     [x,y,b]=ginput(1);
%     x=max(1,min(nim,round(x)));
%     disp([num2str(x) fname{x}]);
% end;
