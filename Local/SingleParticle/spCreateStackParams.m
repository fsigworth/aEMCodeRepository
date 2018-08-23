function sp=spCreateStackParams(stack,pixA)
% Create the sp structure to go with a particle-image stack
[n, ny, nim]=size(stack);
sp=struct;
sp.pixA=pixA;
sp.boxSize=n;
sp.trans=zeros(nim,2,'single');  % in original pixels
sp.rot=zeros(nim,1,'single');    % in degrees
sp.class=zeros(nim,1,'single');
sp.thetaPhi=zeros(nim,2,'single');
sp.cc=zeros(nim,1,'single');
sp.amp=zeros(nim,1,'single');
sp.active=true(nim,1);
sp.flip=false(nim,1);

