% rtCopyJpegTree
% in a directory tree, create a Jpegs folder, then copy into it all jpeg files
% in the tree.  This lets us pick up images from EPU.
% eg given pwd = experiment/Kvdata/ , create Jpegs/  Then copy the
% tree and all *.jpg files into this directory, including Jpeg/Kvdata/...*.jpg

jpegDir='../Jpegs/';

p.extensions={'.jpg'};
p.displayOn=0;
fNames=SearchDirectories('./',cell(0),p);
ok=mkdir(jpegDir);
ok
nf=numel(fNames);
disp([num2str(nf) ' files found.']);
for i=1:nf
    name=fNames{i};
    [pa,nm,ex]=fileparts(name);
    CheckAndMakeDir([jpegDir pa],1);
    str=['cp ' name ' ' jpegDir name];
%     disp(str);
    system(str);
end;
