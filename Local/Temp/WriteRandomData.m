% WriteRandomData

for j=1:10
    m=8*rand(4096,4096,8);
    fname=['x' num2str(j) '.tif']
    beep
    WriteTiffStack(m,1,fname);
    disp('done');
    beep; pause(0.2);beep
end;