% StackOutlierEditor.m
% Removes particles from a stack having more than threshSD power in the
% frequency range from 1/160 to 1/80 angstroms.  The particle images and
% entries in the si structure are deleted outright.  The modified stack
% and si structure are written out with the outPrefix added to the name.

outPrefix='e1';
d0=160;  % angstroms min freqency
d1=80; % angstroms
threshSD=2.5;

[fname, stackPath]=uigetfile('*.mat','Select tsi files','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
cd(stackPath);

nst=numel(fname);
%%
for index=1:nst
stName=fname{index};
    q1=strfind(stName,'tsi.mat');
    if numel(q1)>0 && q1(end)>1
        load(stName);
        si0=si;  % safe copy
        baseName=stName(1:q1(end)-1);
        stackName=[stackPath baseName 'stall.mrc'];
        disp(['Reading ' stackName]);
        stack=ReadMRC(stackName);
        
        [nx ny nim]=size(stack);
        numImages=nim
        sps=RadialPowerSpectrum(stack);
        fmin=floor(si.pixA*nx/d0);
        fmax=ceil(si.pixA*nx/d1);
        lfv=(fmin:fmax)*sps(fmin:fmax,:);
        sigma=std(lfv);
        maxVar=threshSD*sigma
        figure(3);
        hist(lfv,100);
        xlabel('LF variance');
        ylabel('Number of particles');
        excl=lfv>threshSD*sigma;
        numExcluded=sum(excl)
        
        stackEd=stack;
        stackEd(:,:,excl)=0;
        minStackEd=min(stackEd(:));
%         stackEd(:,:,excl)=minStackEd;

        figure(5);
        ImagicDisplay2(BinImage(stack,2));
        figure(6);
        ImagicDisplay2(BinImage(stackEd,2)); %%
        
%%         remove particles from stack file and info file
        
stackShort=stack;
        stackShort(:,:,excl)=[];
        
        fields=fieldnames(si);
        nmod=0;
        for i=1:numel(fields)
            q=si.(fields{i});
            if numel(q)==nim  % this array counts each image
                q(excl)=[];
                si.(fields{i})=q;
                nmod=nmod+1;
            end;
        end;
        disp([num2str(nmod) ' si fields updated.']);
        
        if nmod>0
            outName=[stackPath baseName outPrefix];
            disp(['Writing ' outName 'tsi.mat']);
            save([outName 'tsi.mat'],'si');
            disp(['Writing ' outName 'stall.mrc']);
            WriteMRC(stackShort,si.pixA,[outName 'stall.mrc']);
        end;
        
    end;
end;
%% 

figure(4);
SetGrayscale;
imacs(sum(stackEd,3));
