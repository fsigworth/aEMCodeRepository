function [spec, shot, ok]=meEvalNoiseModel(f,mi)
% function [spec shot]=meEvalNoiseModel(f,mi)
% Given the noise model information stored in the micrograph info
% structure, compute the specimen spectrum (spec) and the CTF-independent
% spectrum (shot) at the given frequencies f (A^-1)
% We evaluate the code stored in mi.noiseModelCode.  This 
% code expects the input variables f and p, and returns spec and shot.
spec=zeros(size(f));
shot=zeros(size(f));
ok=true;
p=[];
if isfield(mi,'noiseModelCode') && isfield(mi,'noiseModelPars')
    nc=numel(mi.noiseModelCode);
    p=mi.noiseModelPars;
end;
if numel(p)<1 || (nc<1)
    if nargout<3
        error('Invalid noise model');
    else
        ok=false;
        return;
    end;
end;
for i=1:nc
    eval(mi.noiseModelCode{i});
end;
