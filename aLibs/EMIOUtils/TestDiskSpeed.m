% TestDiskSpeed
% Try writing and reading a big file

[nm,pa]=uiputfile('*.x','File to write and read','test.x');
cd(pa);
f=fopen(nm,'w');

dato=rand(8192,8192,2);  % 1 GB of data

disp('Writing');
tic
fwrite(f,dato,'double');
toc
fclose(f);

f=fopen(nm,'r');
disp('Reading');
tic
dati=fread(f,'*double');
toc
fclose(f);
