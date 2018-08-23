function mi=meStoreNoiseModel(p,functionName,mi)
% function mi=meStoreNoiseModel(p,functionName,mi)
% Copy the text of the m-file for <functionName>
% into a cell array in the mi structure; also
% copy the parameter vector p into the structure.  The noise model function
% can then be evaluated by meEvalNoiseModel().
% We delete the first line of the m-file text to avoid
% the function declaration. 
% The new mi fields are
% mi.noiseModelCode  (cell array of strings)
% mi.noiseModelPars  (vector of floats or doubles)

fullName=which(functionName);
if ~FileExists(fullName)
    error(['Can''t find the file ' fullName]);
end;
f=fopen(fullName);
q=textscan(f,'%s','delimiter','');
fclose(f);
qc=q{1};  % pick up the first element, which is a cell array
nc=numel(qc);
if nc<2
    error(['Not a valid m-file, too short: ' fullName]);
end;
mi.noiseModelCode=qc(2:nc); % Skip the first line.
mi.noiseModelPars=p;
