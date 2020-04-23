function [name,ok]=CheckForAltImage(name,sufExts)
% Check to see if variants of an image with the given name exists.  Construct
% names where sufExts{i} are added to the base of name.  Return the first
% version that is found to exist.
% For example,
%  sufExts={'s.mrc' 'z.tif' '.mrc'}
%  [name,ok]=CheckForAltImage('imagem.mrc',sufExts);
% ... could return for example name='imagems.mrc' if that exists.

    i=0;
    ok=false;
    %     Check to see if there is a compressed version
    [pa,nm,~]=fileparts(name);
    while ~ok && i<numel(sufExts)
        i=i+1;
        name=[AddSlash(pa) nm sufExts{i}];
        ok=exist(name,'file');
    end;
    if ~ok
        if nargout<2
            error(['No file or alternate found with the name ' name]);
        end;
    end;
