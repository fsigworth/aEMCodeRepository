% rlSplitMicrographs


starName='MotionCorr/job002/corrected_micrographs.star';
[nms,dat]=ReadStarFile(starName);
d=dat{1};
names=d.rlnMicrographName(2:end);
nm=numel(names);
outDir='Split816/';
CheckAndMakeDir(outDir);
ns=5; % number of splits
n1=816;
n0=4096;
dx=(n0-n1)/(ns-1)
for i=1:nm
    [m,s]=ReadMRC(names{i});
    [pa,nm,ex]=fileparts(names{i});
    for j=1:ns
          isx=min(n0-n1+1,(j-1)*dx);
        for k=1:ns
            isy=min(n0-n1+1,(k-1)*dx);
            m1=m(isx+1:isx+n1,isy+1:isy+n1);
            outName=[outDir nm '_' num2str(j) num2str(k) ex];
            disp(outName);
            WriteMRC(m1,s.pixA,outName);
        end;
    end;
end;
