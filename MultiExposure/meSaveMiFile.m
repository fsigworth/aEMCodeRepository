function infoName=meSaveMiFile(mi)
% function fullName=meSaveMiFile(mi)
infoName=[mi.infoPath mi.baseFilename 'mi.txt'];
% disp(['saving:  ' fname]);
WriteMiText(mi,infoName);
% save(fullName,'mi');
