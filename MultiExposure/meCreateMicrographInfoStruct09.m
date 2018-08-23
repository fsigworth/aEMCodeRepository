function mi=meCreateMicrographInfoStruct

% Create the MicrographInfo structure 
mi.version=09;
mi.baseFilename='';
mi.pixA=0;
mi.nPix=int16(0);
mi.boxSize=int16(0);
mi.nExposures=int16(0);
mi.camera='';
mi.quality=uint8(zeros(1,8));
mi.notes='';
mi.localPath='';
mi.rawFile={};
mi.doses=0;
mi.ctf=struct([]);
mi.mergeMatrix=[];
mi.vesicleModel=[];
mi.vesicle.x=[];
mi.vesicle.y=[];
mi.vesicle.r=[];
mi.vesicle.amp=[];
mi.particle.x=[];
mi.particle.y=[];
mi.particle.vesicle=[];
mi.particle.quality=[];


% Example of typical values:
% mi.version	10
% mi.baseFilename	?10dec22a_bkoWA451D_00015gr_00042sq_v01_00002hl_v02_00002?
% mi.pixA	1.2			% or should we use pixnm?
% mi.nPix	4096
% mi.boxSize	80
% mi.nExposures	3
% mi.camera	?f20?	% Allows choice of correction parameters, prewhitening
% mi.quality	[0 0 0 0 0 0 0 0]  % eight uint8 values; 0 means unclassified
% mi.notes	?Good ice, but much contamination?
% mi.localPath	?10jan19b17g/?
% mi.rawFile(i)	i=1..3; first element could be ?10dec22a_bkoWA451D_00015gr_00042sq_v01_00002hl_v02_00002ef.mrc?
% mi.doses      [10 20 20] e/A^2
% mi.ctf(i)		% i=1..3
% mi.mergeMatrix(i)		% merge-transformation matrices
% mi.vesicle.model		% 1D scattering profile, in units of rad/A
% mi.vesicle.x(j)		% vector of vesicle x-coords, original pixel units
% mi.vesicle.y(j)
% mi.vesicle.r(j)
% mi.vesicle.amp(j)		% vector of model scalings
% mi.particle.x(l)		% x coord (orig pixels) of particle l in session k
% mi.particle.y(l)
% mi.particle.vesicle(l)	% index (j) of the vesicle containing this particle
% mi.particle.quality(l)	% quality of particle: eight uint8 values.
% 
