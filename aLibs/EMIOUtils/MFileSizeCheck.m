% MFileSizeCheck.m

% First, do a recursive search for *.m files.
% pars.displayOn=0;
% pars.extensions={'.m'};
% pars.maxDepth=4;
% pars.regexps={};

% fnames=SearchDirectories('~/aEMCodeRepository/',{},pars);
fnames=SearchDirectories('');  % take all the defaults
% Now do the correction, listing those files which were modified.
% doWrite=1;
% for i=1:numel(fnames),name=fnames{i}; [conv,ncr,ncrlf]=FixCRsInMFiles(name,doWrite);
%     if conv disp([num2str([ncr ncrlf]) ' ' name]); end;
% end;
