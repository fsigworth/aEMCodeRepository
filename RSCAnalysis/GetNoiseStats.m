% GetNoiseStats

nm=numel(mis);
f=[0 .01];  % have to have at least two elements.
nf=numel(f);
spec=zeros(nm,nf);
def=zeros(nm,1);
datNames=cell(nm,1);
for i=1:nm
    mi=mis{i};
    datNames{i}=[mi.baseFilename 'ms.mrc'];
    if numel(mi.noiseModelPars)>0
        sp=meEvalNoiseModel(f,mis{i});
        spec(i,:)=sp;
        def(i)=mi.ctf(1).defocus;
    else
        spec(i,:)=NaN;
        def(i)=NaN;
    end;
end;
%% bin according to holes
spm=max(spec(:,2),.1);
subplot(2,2,1);
semilogy([spm def/10]);

nh=floor(nm/7);
nm0=nh*7;
spm=spm(1:nm0);
def=def(1:nm0);
spa=reshape(repmat(mean(reshape(spm,7,nh),1),7,1),nm0,1);
dfa=reshape(repmat(mean(reshape(def,7,nh),1),7,1),nm0,1);
plot([spa(:) spm dfa/20 def/20])

return
%% Prepare for selecting images
dfOk=dfa(1,:)'<5;
dfb=dfa(1,dfOk)';  % defocus of our selected set
spb=spa(1,dfOk)';  % spectrum of our selected set


hGood=false(1,nh);
spok=spm;
i=1;  % hole pointer
%% Interactive selection
i=1;
noImageDisplay=1;
fileIndex=0;
if noImageDisplay
    ndis=100*7; % half-display size
else
    ndis=40*7;
end;
% load hGood
b='z';
figure(1);
clf;
subplot(221);
while b~='q'
    i7=7*(i-1)+4;  % micrograph pointer
    rl=max(i7-ndis,4);
    ru=min(rl+2*ndis-1,nm0-3);
    rl=max(ru-2*ndis,4);
    range=rl:ru;
if noImageDisplay
    subplot(1,1,1);
else
    subplot(2,2,1);
end;
    semilogy(range,[spm(range) spa(range)],'-');
    axis([range(1) range(end) .1 3]);
    hold on;
    semilogy(range,spok(range),'r-');
    plot(i7,spa(i7),'ko','markersize',12);
    hold off;
if ~noImageDisplay
    subplot(2,2,2);
    name=['Merged/' datNames{i7}];
    if exist(name,'file')
        imags(ReadMRC(name));
    else
        cla;
    end;
    title(name,'interpreter','none');
    for j=i7-3:i7+3
        if j~=i7
            subplot(4,4,j-i7+12)
            name=['Merged/' datNames{j}];
            if exist(name,'file')
                imags(ReadMRC(name));
                axis off;
            else
                cla;
            end;
        end;
    end;
end;    
    
    [x,y,b]=Myginput(1,'arrow');
    b=char(b);
    switch b
        case char(1)
            i=round((x-4)/7)+1;
        case 'y'
            hGood(i)=true;
            i=min(nh,i+1);
        case 'n'
            hGood(i)=false;
            i=min(nh,i+1);
        case ' '
            i=min(nh,i+1);            
        case 'b'
            i=max(1,i-1);
        case 's'  % save the figure
            fileIndex=fileIndex+1;
            fname=['noise' num2str(fileIndex)]
            print('-r300',fname,'-djpeg');
        otherwise
            beep;
    end;
    good=reshape(repmat(hGood,7,1),nm0,1);
    spok(good)=spm(good);
    spok(~good)=NaN;
    
end;
ok=input('Save the hGood file [1 or 0]? ');
if ok
save hGood hGood dfOk i
disp('hGood.mat saved.');
end;
disp('done.');
return

%%
goodHoles=hGood' & dfOk;
goodImgs=reshape(repmat(goodHoles',7,1),nm0,1);

goodMis=mis(goodImgs);
CheckAndMakeDir('Info_good');
for i=1:numel(goodMis);
        mi=goodMis{i};
        infoName=['Info_good/' mi.baseFilename 'mi.txt'];
        WriteMiText(mi,infoName);
        disp(infoName);
end;


