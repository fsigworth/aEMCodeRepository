function si=rsCreateStackInfoStruct11(nim,nmic)

si=struct;

si.dataFilenames=cell(0,1);
si.urParticle=  zeros(nim,1,'uint32');
si.origParticle=zeros(nim,1,'uint32');
si.miIndex=     zeros(nim,1,'uint16');
si.miParticle=  zeros(nim,1,'uint16');
si.alpha0=      zeros(nim,1,'single');
si.yClick=      zeros(nim,1,'single');  % in units of si.pixA
si.rVesicle=    zeros(nim,1,'single');
si.sVesicle=    zeros(nim,1,'single');
si.mi=          cell(nmic,1);
si.ctfs=        zeros(0,0,nmic,'single');

si.activeFlags= true(nim,1);
si.activeFlagLog={date};
