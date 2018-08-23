% ZTiffReadInfo.m

ok=1;
while ok
% Have the user select a z.tif file
    [fname, pa]=uigetfile('*z.tif','Select a ZTiff file');
    if isnumeric(fname)
        return
    end;
    cd(pa);
    filename=fname;
    [m,s,ok]=ReadZTiff(filename,1,0);
    disp(filename);
    s
    ok=(input('Load another file [n]','s')=='y');
end;
