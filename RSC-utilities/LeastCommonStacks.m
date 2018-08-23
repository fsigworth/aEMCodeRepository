% LeastCommonStacks.m
% Find the particles common to two stacks
epsi=4;
doWrite=0;

% cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Stack');
cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/StackMerge/')

stNames{1}='sq10_350p128tsi.mat';
% stNames{1}='sq10_350nosubtsi.mat';
stNames{2}='sq10_350fp128tsi.mat';

si1=load(stNames{1});
si1=si1.si;
si2=load(stNames{2});
si2=si2.si;

disp('stCommonParticles');

disp('---start----')
disp(stNames{1});
disp(' active flags:');
naf1=size(si1.activeFlags,2); 
for i=1:naf1
    if size(si1.activeFlagLog,1)<i
        si1.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si1.activeFlags(:,i))]) '  ' si1.activeFlagLog{i,1}]);
end;
afi1=size(si1.activeFlags,2);  % Active flag index for si1
afi1=MyInput('flag ',afi1);
active1=si1.activeFlags(:,afi1);

disp(' ');

disp(stNames{2});
disp(' active flags:');
naf2=size(si2.activeFlags,2); 
for i=1:naf2
    if size(si2.activeFlagLog,1)<i
        si2.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si2.activeFlags(:,i))]) '  ' si2.activeFlagLog{i,1}]);
end;
afi2=size(si2.activeFlags,2);
 afi2=MyInput('flag ', afi2);
active2=si2.activeFlags(:,afi2);
disp(' ');



% pick up the corresponding image filenames

np=size(si1.miIndex,1);
nmi1=numel(si1.mi);
mcNames1=cell(nmi1,1);
for i=1:nmi1
    mi=si1.mi{i};
    if numel(mi)>0 && isfield(mi,'imageFilenames')
        mcNames1{i}=mi.imageFilenames{1};
    else
        mcNames1{i}='';
    end;
end;

nmi2=numel(si2.mi);
mcNames2=cell(nmi2,1);
for i=1:nmi2
    mi=si2.mi{i};
    if numel(mi)>0 && isfield(mi,'imageFilenames')
        mcNames2{i}=mi.imageFilenames{1};
    else
        mcNames2{i}='';
    end;
end;

% Make the mcMatch matrix showing where the filenames match
mcMatch=false(nmi1,nmi2);
for j=1:nmi2
    mcMatch(:,j)=strcmp(mcNames1,mcNames2{j});
end;

% Scan through si1, finding all particles that exist in si2 as well.

np1=numel(si1.miIndex);
pIndex2=zeros(np1,1);  % index of corresponding particle in si2

for i=find(active1)' % only search the active particles
    %     First, match up the micrographs
    miInd1=si1.miIndex(i);
    miInd2=find(mcMatch(miInd1,:),1); % Find the corresponding si2.mi index
    if miInd1>0 && numel(miInd2)>0  % we've found a matching micrograph
        mi1=si1.mi{miInd1};
        mi2=si2.mi{miInd2};
%         if isfield(mi1,'particle') && isfield(mi1.particle,'picks') && size(mi1.particle.picks,1)>0
        %         Look for a particle with the same coordinates
 if size(mi1.particle.picks,1)<si1.miParticle(i) || size(mi1.particle.picks,2)<2 ...
         || size(mi2.particle.picks,1)<1
     disp('bounds')
     break
 end;
        xy=mi1.particle.picks(si1.miParticle(i),1:2);               % get the mi1 particle coords
        dists=sqrt((xy(1,1)-mi2.particle.picks(:,1)).^2 + abs(xy(1,2)-mi2.particle.picks(:,2)).^2);
        %  but need correct flags too
        flagsOk2=mi2.particle.picks(:,3)>=16 & mi2.particle.picks(:,3)<48;
        q=find(dists<epsi & flagsOk2);  % mi particle indices tht match
%         else
%             q=[];
%         end;
        %         Check the active flags for these particles in si2
        if numel(q)>1
            ac2s=false(numel(q));
            for j=1:numel(q)
                p=find(si2.miIndex==miInd2 & si2.miParticle==q(j));
                if numel(p)>0
                    ac2s(j)=active2(p);
                end;
            end;
            q=find(ac2s);
            %              warning(['Duplicate particles. si1 particle ' num2str(i)...
            %                 ' flag ' num2str(mi1.particle.picks(si1.miParticle(i),3))...
            %                 ' aflag ' num2str(si1.activeFlags(i,end))...
            %                 ' matches si2 micrograph ' num2str(miInd2) ' particles ' num2str(q')]);
            %             i
        end;
        if numel(q)>0
            p=find(si2.miIndex==miInd2 & si2.miParticle==q(1));
            if numel(p)>0
            pIndex2(i)=p(1);
            end;
        end;
        
    end;
end;

%%

disp('-----end-----');

% Coordinate the active flags
% Create new active flags sets that correspond
mappedAF1=false(np1,1);
inds=pIndex2(pIndex2>0);
mappedAF1(pIndex2>0)=si2.activeFlags(inds,afi2);
mappedAF1=mappedAF1 & si1.activeFlags(:,afi1);
% sum(mappedAF1)
naf1=size(si1.activeFlags,2);
si1.activeFlags(:,naf1+1)=mappedAF1;
si1.activeFlagLog{naf1+1,1}=['stCommonParticles ' date ' with ' stNames{2}];
disp(['Active flags for ' stNames{1}]);
for i=1:naf1+1
    disp([num2str([i sum(si1.activeFlags(:,i))]) '  ' si1.activeFlagLog{i,1}]);
end;
changed1=any(si1.activeFlags(:,afi1)~=si1.activeFlags(:,naf1+1));
disp(' ');


np2=size(si2.miIndex,1);
mappedAF2=false(np2,1);
for i=1:np2
    q=find(pIndex2==i);
    if numel(q)>1
         warning([num2str([i q'])]);
    elseif numel(q)>0
        mappedAF2(i)=si1.activeFlags(q,afi1);
    end;
end;
mappedAF2=mappedAF2 & si2.activeFlags(:,afi2);
% sum(mappedAF2)
disp(['Active flags for ' stNames{2}]);
naf2=size(si2.activeFlags,2);
si2.activeFlags(:,naf2+1)=mappedAF2;
si2.activeFlagLog{naf2+1,1}=['stCommonParticles ' date ' with ' stNames{1}];
for i=1:naf2+1
    disp([num2str([i sum(si2.activeFlags(:,i))]) '  ' si2.activeFlagLog{i,1}]);
end;
changed2=any(si2.activeFlags(:,afi2)~=si2.activeFlags(:,afi2+1));
disp(' ');

if doWrite && (changed1 || changed2)
    disp('Writing new files');
    si=si1;
    save(stNames{1},'si');
    si=si2;
    save(stNames{2},'si');
else
    disp('Nothing written');
end;

