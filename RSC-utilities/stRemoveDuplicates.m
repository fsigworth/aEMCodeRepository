function stRemoveDuplicates

[siName, pa]=uigetfile('*si.mat','Select si file');
if isnumeric(pa)
    return
end;
cd(pa);
disp('Reading:');
disp(siName);
si=LoadStruct(siName);

disp(' ');
disp(' active flags:');
naf1=size(si.activeFlags,2);
for i=1:naf1
    if size(si.activeFlagLog,1)<i
        si.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si.activeFlags(:,i))]) '  ' si.activeFlagLog{i,1}]);
end;
k=naf1;
k=MyInput('Select data at selection ',k);
active=si.activeFlags(:,k);
initialActive=sum(active)
disp('Deleting particles:');
np=size(si.miIndex,1);
mInd=0;
inds={};
flags={};
for ip=1:np
    if si.miIndex(ip)~=mInd
        mInd=si.miIndex(ip);
        [indSet,flagSet]=rspFindDuplicateParticles(si.mi{mInd});
        if numel(indSet)>0 % we have duplicates
            for di=1:numel(indSet)
                inds=indSet{di}; % indices of a set in particles in this micrograph
                flags=flagSet{di};

                pinds=ParticleLookup(si,mInd,inds);
                pOk=pinds>0;
                pindso=pinds(pOk);
                flags=flags(pOk);
                pActive=active(pindso);
                if sum(pActive)>1 % we have duplicates
                    pindsa=pindso(pActive);  % particles that are active and duplicates
                    flags=flags(pActive);
                    [fv,finds]=sort(flags);
%                 We'll keep the highest flag value but delete the rest.
                    finds(end)=[];
                    disp(pindsa(finds));
                    active(pindsa(finds))=0;
                end;
            end;
        end;
    end;
end;
disp('done.');
finalActive=sum(active)
if finalActive<initialActive
    naf=size(si.activeFlags,2);
    si.activeFlags(:,naf+1)=active;
    si.activeFlagLog{naf+1,1}=['stRemoveDuplicates (' num2str(k) ') ' date];
    disp(['Writing ' siName]);
    save(siName,'si');
else
    disp('No change, nothing written');
end;
                
                
                function partInds=ParticleLookup(si,miIndex,miParticles)
                nps=numel(miParticles);
                partInds=zeros(1,nps);
                for i=1:nps
                    q=find(si.miIndex==miIndex & si.miParticle==miParticles(i));
                    if numel(q)>0
                        partInds(i)=q(1);
                    end;
                end;
                