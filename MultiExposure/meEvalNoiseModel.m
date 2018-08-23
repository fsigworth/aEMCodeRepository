function [spec, shot]=meEvalNoiseModel(f,mi)
% function [spec shot]=meEvalNoiseModel(f,mi)
% Given the noise model information stored in the micrograph info
% structure, compute the specimen spectrum (spec) and the CTF-independent
% spectrum (shot) at the given frequencies f (A^-1)
% We evaluate the code stored in mi.noiseModelCode.  This 
% code expects the input variables f and p, and returns spec and shot.
nc=numel(mi.noiseModelCode);
p=mi.noiseModelPars;
if numel(p)<1 || (nc<1)
    error('Invalid noise model');
end;
for i=1:nc
    eval(mi.noiseModelCode{i});
end;
