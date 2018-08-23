% rsUpdateKvRsos
% Update particle picks to include the rso information.  This is for
% particles picked before 23 July 2014.
doWrite=1;  % Actually write out the new mi file.
% Have the user select some mi files: boilerplate
[fname, pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
if isnumeric(pa) % File selection cancelled
    return
end;
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);
%%
for fileIndex=1:numel(fname); % Operate on a single micrograph
    miName=[infoPath fname{fileIndex}];
    disp(['Reading ' miName]);
    mi=ReadMiFile(miName);
    mi.basePath=rootPath;
    %     See if we have anything to do
    if isfield(mi,'particle') && isfield(mi.particle,'picks')...
            && size(mi.particle.picks,1)>0
        picks=mi.particle.picks;
        nPicks=size(picks,1);
        rsccName=[mi.procPath mi.baseFilename 'rscc.mat'];
        if exist(rsccName,'file')
            rscc=load(rsccName);
            dsc=mi.imageSize(1)/size(rscc.mxCC,1);  % Downsample factor of mxCC
            oldRsos=0;
            newRsos=0;
            for i=1:nPicks
                if picks(i,3)>=16 % a picked particle
                    dsCoords=round(picks(i,1:2)/dsc)+1;  % downsampled coordinates to look up attributes
                    oldRsos=oldRsos+picks(i,7);
                    picks(i,7)=single(rscc.mxRsos(dsCoords(1),dsCoords(2)));
                    newRsos=newRsos+picks(i,7);
                end;
            end;
            disp([rsccName '   old, new Rsos: ' num2str([oldRsos newRsos])]);
            if doWrite
                mi.particle.picks=picks;
                if ~isfield(mi,'log')
                    mi.log=[];
                end;
                mi.log{end+1,1}=['rsUpdateKvRsos ' TimeStamp];
                WriteMiFile(mi,miName);
            end;
        else
            warning(['Couldn''t find ' rsccName]);
        end;
    else
        disp('no picks.');
    end;
end;