function pars=SetOptionValues(defaults,newPars,warningsOn)
% function pars=SetOptionValues(defaults,newPars,warningsOn)
% Make the struct pars with fields from defaults except when there are 
% corresponding fields in newPars, which then are used.  We call this giving the
% defaults and user-provided newPars, to get the final parameters for a function.
% By default, warnings are give in a field in newPars is not present in defaults.

if nargin<3
    warningsOn=1;
end;

pars=defaults;
parFields=fieldnames(pars);
newFields=fieldnames(newPars);
for i=1:numel(newFields)
    p=strcmp(newFields{i},parFields);
    if any(p)
        ptr=find(p,1);
        pars.(parFields{ptr})=newPars.(newFields{i});
    elseif warningsOn
        disp(['  Unrecognized option: ' newFields{i}]);
    end;
end;
