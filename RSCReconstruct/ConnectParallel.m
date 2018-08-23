jmName = fgetl(fopen('MDCE_JM_NAME')) 
jm = findResource('scheduler', 'type', 'jobmanager', 'Name', jmName)
matlabpool(jm, 'open')
matlabpool('size')