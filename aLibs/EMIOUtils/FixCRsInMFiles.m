function [converted,ncr,ncrlf]=FixCRsInMFiles(filename,doWrite);
% Find CR characters (from old MacOS) and convert them to LF characters
% (Unix standard).  This avoids a bug in Matlab where the help feature
% doesn't recognize CRs as newline characters.

% %you can fix all the m files thusly:
%
% % First, do a recursive search for *.m files.
% pars.displayOn=0;
% pars.extensions={'.m'};
% pars.maxDepth=4;
% fnames=SearchDirectories('~/aEMCodeRepository/',{},pars);
% % Now do the correction, listing those files which were modified.
% doWrite=1;
% for i=1:numel(fnames),name=fnames{i}; [conv,ncr,ncrlf]=FixCRsInMFiles(name,doWrite);
%     if conv disp([num2str([ncr ncrlf]) ' ' name]); end;
% end;


h=fopen(filename);
s=fread(h,'*uint8');
fclose(h);
ns=numel(s);
p=find(s==13);  % find CR characters
ncr=0;
ncrlf=0;
if numel(p)>0
    for i=1:numel(p)
        ptr=p(i);
        if ptr<ns && s(ptr+1)==10 % we have a CRLF: don't do anything
            ncrlf=ncrlf+1;
            continue;
        else
            s(ptr)=10;  % convert to LF
            ncr=ncr+1;
        end;
    end;
end;
converted=(ncr>0);
if converted && doWrite
    ho=fopen(filename,'w');
    fwrite(ho,s);
    fclose(ho);
end;

