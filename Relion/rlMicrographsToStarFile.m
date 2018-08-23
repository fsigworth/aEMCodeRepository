% rlSiToStarFile
% Create a .star file from our mi file, pointing to the 1st exposure
% micrographs.
% 
% data_images
% loop_
% _rlnImageName
% _rlnDefocusU
% _rlnDefocusV
% _rlnDefocusAngle
% _rlnVoltage
% _rlnAmplitudeContrast
% _rlnSphericalAberration
% ...
% 000001@/sq02_1n96tstack.mrc 13538 13985 109.45 300 0.15 2
% 000002@/sq02_1n96tstack.mrc 13293 13796 109.45 300 0.15 2
% 000003@/sq02_1n96tstack.mrc 13626 14085 109.45 300 0.15 2
dataName='micrographs';
detPixel=5;
ampContrast=.1;

% Put up a file selector for info files
disp('Getting mi files.');
[miName, miPath]=uigetfile('*mi.txt','Select mi files','multiselect','on');
if isnumeric(miPath)
    return
end;
if ~iscell(miName)
    miName={miName};
end;
nmi=numel(miName);
disp([num2str(nmi) ' mi files.']);
%%
% Output file name
[starName, starDir]=uiputfile('.star','Pick an output star file');
if ~ischar(starDir)
return
end;

cd(starDir);
disp(['Output file: ' starName]);
%


fi=fopen(starName,'w');
% fi=1;

fprintf(fi,'data_%s\n',dataName);

fprintf(fi,'loop_\n');
fprintf(fi,'_rlnMicrographName\n');
fprintf(fi,'_rlnDefocusU\n');
fprintf(fi,'_rlnDefocusV\n');
fprintf(fi,'_rlnDefocusAngle\n');
fprintf(fi,'_rlnVoltage\n');
fprintf(fi,'_rlnAmplitudeContrast\n');
fprintf(fi,'_rlnSphericalAberration\n');
fprintf(fi,'_rlnMagnification\n');
fprintf(fi,'_rlnDetectorPixelSize\n');

%allActive=activeFlags & sdOk;
for i=1:nmi
        mi=ReadMiFile([AddSlash(miPath) miName{i}]);
        imgName=[mi.imagePath mi.imageFilenames{1}];
        mag=detPixel*1e4/mi.pixA;
%        mName=[mi.baseFilename 'm.mrc'];
        ctf=mi.ctf(1);
        defocus=ctf.defocus;
        deltadef=ctf.deltadef;        
        defU=round((defocus+deltadef)*1e4);
        defV=round((defocus-deltadef)*1e4);
        ang=mod(ctf.theta*180/pi,180);
        kV=mi.kV;
        Cs=ctf.Cs;
            fprintf(fi,'%s %d %d %6.2f %g %5.3f %g %g %g\n',imgName,defU,defV,ang,kV,ampContrast,Cs,mag,detPixel);
end;
fprintf(fi,'\n');
fclose(fi);
disp(' written.');