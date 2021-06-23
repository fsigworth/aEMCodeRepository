function C=StructsConcat(A,B)
% function C=StructsConcat(A,B)
% For every field common to A and B, e.g. A.nm and B.nm, return the 
% concatenated column vectors C.nm=[A.nm ; B.nm];
% 
% fields are assumed to be column vectors.
fa=fieldnames(A);
fb=fieldnames(B);
C=struct;
na=numel(fa);
for i=1:na
    fname=fa{i};
    if isfield(B,fname)
        C.(fname)=[A.(fname) ; B.(fname)];
    end;
end;

