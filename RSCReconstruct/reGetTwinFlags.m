function flags=reGetTwinFlags(ri,nActive,iTwin)
% Get the binary flags for the twin collection of images.
% nActive is the number of images in the already squeezed stack.
% it iTwin is given as zero, flags is returned as an nActive x 2 array of
% booleans giving flags for both twins.

% Figure out the twin flags
if iTwin==0
    nt=ri.nTwins;
else
    nt=1;
end;
flags=false(nActive,nt);
for it=iTwin:iTwin+nt-1
    flags(it:ri.nTwins:end)=true;
end;
