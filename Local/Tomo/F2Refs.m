% F2Refs.m


cd('/Users/fred/EMWork/Hideki/Falcon2ref')
d=dir;

m=ReadEMFile(d(13).name);
for j=2:5
    m(:,:,j)=ReadEMFile(d(j+12).name);
end;
mf=single(m);


