function reSaveMat(fileName,var,varName)

tmpName=[fileName '_tmp'];
save(tmpName,var);
pause(1);
cmdString=['!mv ' tmpName ' ' fileName'];
exec(cmdString);
