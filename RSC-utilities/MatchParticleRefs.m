% MatchParticleRefs.m
% Given two si files si1 and si2, compare the picked particles by micrograph and location.
% Create a new stack with images from si2 but in the order of si1.  Missing images
% are indicated by the blanks boolean vector and are filled with noise.
% I use this to coordinate single-exposure images with merged image stack.

[file1, pa1]=uigetfile('*si.mat','Select the reference si file');
if isnumeric(pa1)
    return
end;
cd(ParsePath(pa1));  % back up by one directory
pa1=AddSlash(pa1);

[file2, pa2]=uigetfile('*si.mat','Select si file to modify');
if isnumeric(pa2)
    return
end;
pa2=AddSlash(pa2);
cd(pa2);

disp(['Reading ' pa1 file1]);
si1=load([pa1 file1]);
si1=si1.si;
nmi1=numel(si1.mi);
mid1=cell(nmi1,1);
for i=1:nmi1
    mid1{i}=si1.mi{i}.baseFilename;
end;

disp(['Reading ' pa2 file2]);
si2=load([pa2 file2]);
si2=si2.si;
nmi2=numel(si2.mi);
mid2=cell(nmi2,1);
for i=1:nmi2
    mid2{i}=si2.mi{i}.baseFilename;
end;

%%
% Get the mapping: the si2 micrograph corresponding to si1.
ind2=zeros(nmi1,1);
for i=1:nmi1
    q=find(strcmp(mid1{i},mid2),1);
    if numel(q)>0
        ind2(i)=q;
    end;
end;

disp([num2str(sum(ind2>0))]);
% Find matching particles
% for each particle, find its mi2 file.
% see if a corresponding particle is listed.  Check its location

np=numel(si1.miIndex);
pInd2=zeros(np,1); % pointer in stack2 of ith stack1 particle.
dThresh=1;  % allow 4 pixel error.
count1=0;
count2=0;
for i=1:np
    miInd=si1.miIndex(i);
    if ind2(miInd)>0 % matching mi exists
        mi1=si1.mi{miInd};
        mi2=si2.mi{ind2(miInd)};
        mip=si1.miParticle(i);
        pc1=mi1.particle.picks(mip,1:3);
        pc2s=find(si2.miIndex==ind2(miInd));  % list of stack particles from si2
        if numel(pc2s)>0
            count1=count1+1;
            pm2s=si2.miParticle(pc2s);  % mi indices of the extant particle images
            coords2=mi2.particle.picks(pm2s,1:3);
            dists=sqrt(sum((coords2-repmat(pc1,numel(pc2s),1)).^2,2));
            [minDists,dInd]=min(dists); % get the squared distances
            if minDists<dThresh
                pInd2(i)=pc2s(dInd);
                count2=count2+1;
            end;
        end;
    end;
end;
count1
count2
%%
% Now, read the stack
inSuffix={'stack.mrcs' 'ustack.mrcs'};
% outSuffix={'stack.mrcs_match' 'ustack.mrcs_match'};
outSuffix={'stack_match.mat' 'ustack_match.mat'};
baseName=file2(1:end-6);
blanks=pInd2==0;

for k=1:2
    stack2Name=[baseName inSuffix{k}];
    if exist([pa2 stack2Name],'file')
        disp(['Reading ' stack2Name]);
        [m2,s]=ReadMRC([pa2 stack2Name]);
        n=size(m2,1);
        disp('Sorting the stack');
        m2New=zeros(n,n,np,'single');
        for i=1:np
            if pInd2(i)>0
                m2New(:,:,i)=m2(:,:,pInd2(i));
            end;
        end;
        
        disp(['Inserting ' num2str(sum(blanks)) ' blanks.']);
        m2New(:,:,blanks)=randn(n,n,sum(blanks));
        outName=[pa2 baseName outSuffix{k}];
        disp(['Saving ' outName]);
        save(outName,'m2New','pInd2','s','-v7.3');
    else
        disp(['Doesn''t exist: ' pa2 stack2Name]);
    end;
end;
clear m2 m2New

%

% mi1Particles=cell(nmi1,1);
% mi1PartInds=cell(nmi1,1);
% mi2Particles=cell(nmi1,1);
% mi2PartInds=cell(nmi1,1);
% match1=false(nmi1,1);
% for i=1:nmi1
%     mi1Particles{i}=find(si1.miIndex==i);  % particle numbers from this mi
%     mi1PartInds{i}=si1.miParticle(mi1Particles{i}); % The corresponding miParticles
%     mi2Particles{i}=find(si2.miIndex==ind2(i));
%     mi2PartInds{i}=si2.miParticle(mi2Particles{i});
%     match1(i)=numel(mi1Particles{i})==numel(mi2Particles{i}) && all(mi1Particles{i}==mi2Particles{i});
% end;
%
% % a valid particle will have
%
% return
%
%
