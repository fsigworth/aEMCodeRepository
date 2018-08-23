% miCopyVesicleData
% Copy the vesicleModel and vesicle fields of one set of mi files to
% another.  For each source file of a given name, look in the target directory
% for a file of the same name.  If a target file isn't found a warning is issued.
% This script is useful in cases where images have been re-merged e.g.
% to change the downsampling factor of the merged images.  In this case the
% vesicle parameters can remain the same and don't need re-fitting.

testMode=1  % Don't write anything if this is set.
doReplaceMembraneModel=1;  % copy the membrane model (mi.vesicleModel)
doReplaceVesicleFits=0;

% Pattern for filenames
rexp='.+mi\.mat';  % some characters followed by mi.mat

[fname, infoPath]=uigetfile('*mi.mat','Select source mi files','multiselect','on');
if isnumeric(fname)  % Cancel
    return
end;
if ~iscell(fname)
    fname={fname};
end;
cd(infoPath);

% Put up a second file selector
infoPath2=uigetdir('../','Select the directory of target files');
if isnumeric(infoPath2)  % Cancel
    return
end;
infoPath2=AddSlash(infoPath2);
%%

for i=1:numel(fname);
    miName=fname{i};
    q=regexp(miName,rexp);
    if exist(miName,'file') && numel(q)>0
        mi2Name=[infoPath2 fname{i}];
        if exist(mi2Name,'file')
            mi1=load(miName);
            mi1=mi1.mi;
            mi2=load(mi2Name);
            mi2=mi2.mi;
            
            if doReplaceMembraneModel
                mi2.vesicleModel=mi1.vesicleModel;
            end;
            if doReplaceVesicleFits
                mi2.vesicle=mi1.vesicle;
            end;
            mi=mi2;  % we'll write this out
            if testMode
                disp(['Would be modified: ' mi2Name]);
            else
                save(mi2Name,'mi');
                disp(['Updated: ' mi2Name]);
            end;
        else
            warning(['Target file doesn''t exist: ' miName]);
        end;
    elseif ~exist(miName,'file')
        warning(['Source file doesn''t exist: ' miName]);
    end;
end;
