% rlCurateParticles

processTheStack=0;

% ds=2;  % downsample to 96/2=48
nCrop=32;
ndis=1000;
comp=.4;
baseName='sq03_5543m2rp96m2';

if ~exist('si','var')
    disp(['Loading ' baseName 'tsi.mat']);
    load([baseName 'tsi.mat']);
end;

rotStackName=[baseName '_rotstack' num2str(nCrop) '.mat'];
if processTheStack
    disp(['Reading ' baseName 'tstack.mrcs']);
    m=ReadMRC([baseName 'tstack.mrcs']);
    n=size(m,1)/ds;
    
    md=Downsample(m,n,1);
    clear m %%%%%%
    nim=size(md,3);
    % nim=5000;
    ndis=1000;
    disp('Inverse filter');
    m1=rlCTFInverseFilter(md(:,:,1:nim),si,ds,comp,.002);
    disp('Rotation');
    m2=Crop(rsRotateImage(m1,si.alpha0(1:nim)),nCrop,1);
    disp('...done.');
    disp(['Saving ' rotStackName]);
    save(rotStackName,'m2');
elseif ~exist('m2','var')
    load(rotStackName);
end;

% display only the active flag set
afIndex=rlPickActiveFlagSet(si);
disFlags=si.activeFlags(:,afIndex);
m3=m2(:,:,disFlags);

nim=size(m3,3);
if ndis>nim
    ndis=nim;
end;
nCrop=size(m3,1);
stackOffset=0;

actName=[baseName '_allActive.mat'];
if exist(actName,'file')
    disp(['Loading ' actName]);
    load(actName);
else
    disp('New allActive array.');
    allActive=disFlags;
end;

active=allActive(disFlags);  % map to displayed set


    disp(['Displaying ' num2str(nim) ' particles.']);

markers=['ys';'cs';'gs';'ms'];
clickMarker=['rx';'kx'];
mk.markerSize=nCrop;

figure(1);
%%
% imagsar(m3(:,:,stackOffset+1:stackOffset+ndis),.0001);
b='r';

while (b~='q') && (b~='Q') && (b~=4) % q = quit; double-click=quit.
    
    switch b
        case {1 2 3 'x'}  % mouse button  1: toggle bad;
            %             key 'x': mark the whole range bad
            %               you select a range by left-clicking the the first element,
            %               then button 2 or x on the last element.
            stackIndex=clickIndex+stackOffset;
            disp([' particle ' num2str(stackIndex)]);
            %             stackIndex=MyInput(' Particle index',stackIndex);
            if stackIndex>0 && stackIndex <=nim
                % -----------handle mouse buttons ---------------
                if b==1  % toggle active
                    if active(stackIndex)
                        mk.marker=clickMarker(1,:); % mark bad
                        active(stackIndex)=false;
                        disp('not active.');
                    else
                        mk.marker=clickMarker(2,:); % mark active
                        active(stackIndex)=true;
                        disp('active.');
                    end;
                    imagsar('mark',clickIndex,mk);
                    prevPartIndex=stackIndex;
                elseif (b==2 || b=='x') && prevPartIndex>0
                    q=active(prevPartIndex); % shift-click: extend state;
                    for p=(prevPartIndex+1:stackIndex)
                        if q
                            mk.marker=clickMarker(2,:); % activate
                            active(p)=true;
                        else
                            mk.marker=clickMarker(1,:); % erase
                            active(p)=false;
                        end;
                        imagsar('mark',p-stackOffset,mk);
                    end;
                    %                 elseif b==3
                    %                     if active(dis.stackIndex)
                    %                         mk.marker='ro';
                    %                     else
                    %                         mk.marker='ko';
                    %                     end;
                    %                     mk.markerSize=0.6*n/dis.stackds;
                    %                     imagsar('mark',clickIndex,mk);
                    %                     active(dis.stackIndex)=~active(dis.stackIndex);
                end;
            end;
        case {'r','n','N'} % shift or refresh
            if b=='n'
                stackOffset=max(0,min(stackOffset+ndis,nim-ndis));
            elseif b=='N'
                stackOffset=max(0,stackOffset-ndis);
            end;
            disp(stackOffset);
            imagsar(m3(:,:,stackOffset+1:stackOffset+ndis),.0001);
            marks=active(stackOffset+1:stackOffset+ndis);
            for i=find(~marks)'
                mk.marker=clickMarker(1,:); % mark inactive
                imagsar('mark',i,mk);
            end;
            set(gcf,'name',['Particles ' num2str(stackOffset+1) '-' num2str(stackOffset+ndis)]);
    end; % switch
    
    oldB=b;  % store the previous key
    pause(.05);
    [clickIndex,coords,b]=imagsar;
    if numel(b)<1
        b=0;
    end;
end; % while

allActive=rlMapActiveFlags(disFlags,active);
disp([num2str(sum(allActive)) ' active of ' num2str(numel(active)) ' displayed, ' ...
    num2str(numel(allActive)) ' total.']);
save(actName,'allActive','stackOffset','afIndex');
disp([actName ' written.']);

