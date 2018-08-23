% StackCombiner.m
% Merge multiple stacks into one.  We assume that pixA and image size are
% the same among the stacks.  A change in image size will give a fatal error.
% The output files are named like the last input file, but will have a '+',
% e.g.
% XXXX+tsi.mat
% XXXX+stall.mrc

stackSuffix='stall.mrc';
% stackSuffix='tstack.mrc';

ds=2;

% Get a tsi file
[fname, pa]=uigetfile('*tsi.mat','Select a tsi file');
cd(pa);
%%
% Load the first one and use it as reference.
disp(['Reading ' fname]);
load(fname);
nmis=numel(si.mi);
disp([num2str(nmis) ' micrographs']);
nims=numel(si.miIndex);  % number of particles
disp([num2str(nims) ' particles']);

p=strfind(fname,'tsi.mat');
stackName=[fname(1:p-1) 'stall.mrc'];
disp(['Reading ' stackName]);
stack=ReadMRC(stackName);
[n ny nstk]=size(stack);
disp(['Image size: ' num2str(n)]);
if nstk ~= nims
    error(['Inconsistent stack length: ' num2str([nstk nims])]);
end;
outStack=stack;

for j=1:nmis  % loop through the micrographs
    mi=si.mi{j};
    imIndices=find(si.miIndex==j); % particles from this micrograph    
    nIndices=numel(imIndices);
    if nIndices>0
        h=ifftshift(meGetCTFInverseFilter(mi,n,ds));
        for i=1:nIndices
            outStack(:,:,imIndices(i))=...
                real(ifftn(fftn(stack(:,:,imIndices(i))).*h));
        end;
    end;
end;
    outBasename=[fname(1:p-1) 'f'];
    stackOutName=[outBasename 'stall.mrc'];
    WriteMRC(stack,si.pixA,stackOutName);
    disp(['Writing ' stackOutName]);
