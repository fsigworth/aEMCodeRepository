% % % FindIOParticles.m
% % % Run this after loading the stacks in reEMReconstruct2 with the proper
% % % fscFlags.
%     We insert a new column into si.activeFlags.
% % 
% % cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon128c')
% % 
% % 
% % % % fscFlag=-1;
% % % load 141106OutputVols1t.mat
% % 
% % % fscFlag=-2;
% % load 141106OutputVols2t.mat
% % 
% % 
% % iter=numel(rois)
% % r=rois{iter};
% % nAI=size(r.pAlphas,1);
% % pRI=[sum(r.pAlphas(1:nAI/2,:))' sum(r.pAlphas(nAI/2+1:nAI,:))'];
% % flagIO=pRI(:,1)<pRI(:,2);
% % fracIO=sum(flagIO)/numel(flagIO)
% % residIO=sum(pRI(flagIO,1))/sum(pRI(:))
% % ioIndices=si.origParticle(flagIO);
% % 
% % 
% % save ioIndicesT2.mat ioIndices
% % %%
% % ioIndices2=ioIndices;
% % %
% % load ioIndicesT1.mat
% % ioIndices1=ioIndices;
% % %%
% % 
siName=[stackPath names{1}];
% % load(siName);  % get si
% % 
% % nf=size(si.activeFlags,2);
% % si.activeFlags(:,nf+1)=false;
% % si.activeFlags(ioIndices1,nf+1)=true;
% % si.activeFlags(ioIndices2,nf+1)=true;
% % si.activeFlags(:,nf+1)=si.activeFlags(:,nf+1) & si.activeFlags(:,nf);
% % 
% % sum(si.activeFlags)
% % 
% % 
% Compare with particle picker results

aInds=find(si.activeFlags(:,nf));
nTrueIO=0;
nFalseIO=0;
badIp=0;
for i=aInds'
    mi=si.mi{si.miIndex(i)};
    np=size(mi.particle.picks,1);
    ip=si.miParticle(i);
    if ip<=np
        nio=mi.particle.picks(si.miParticle(i),7);
        tio=si.activeFlags(i,nf+1);
        if ~nio && tio
            nTrueIO=nTrueIO+1;
        elseif nio && tio
            nFalseIO=nFalseIO+1;
        end;
    else
        badIp=badIp+1;
    end;
end;
nTrueIO
nFalseIO
badIp


%%
si.activeFlagLog{nf+1}=[date '  FindIOParticles'];
save siName si


