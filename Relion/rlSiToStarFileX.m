% rlSiToStarFile
% Create a .star file from our si.mat file, pointing to the existing .mrc
% stack files.  Typical output:
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
dataName='images';
stackSuffixS='stack.mrc';
stackSuffixU='ustack.mrc';
rlSuffix='s';  % add to stack name

unsubToo=1;  % include the unsubtracted stack
checkSD=1;
setWeights=0; % change weights in mi file copy
minSD=.01;  % criterion for including image
detPixel0=5;     % original pixel size in um
ampContrast=.1;  % if nonzero, override ctf(1).alpha value with this.

%ampContrast=.1;  % if nonzero, override ctf(1).alpha value with this.

activeIndex=inf; % use the last active index set.

% Put up a file selector for *si.mat,
disp('Getting an si file.');
[siName, siPath]=uigetfile('*si.mat','Select si file');
if isnumeric(siPath)
    return
end;
cd(siPath);
disp(siName);
si=load(siName);
si=si.si;
nim=numel(si.miIndex);
afIndex=min(activeIndex,size(si.activeFlags,2));
activeFlags=si.activeFlags(:,afIndex);

fprintf('Active flag set %d, %d active out of %d.\n',afIndex,sum(activeFlags),numel(activeFlags));

% Construct the stack data names
[pa,nm,ex]=fileparts(siName);
if numel(nm)<3
    error('si name is too short');
end;
sDatName=[nm(1:end-2) stackSuffixS];
uDatName=[nm(1:end-2) stackSuffixU];
sDatNameRl=[nm(1:end-2) stackSuffixS rlSuffix];
uDatNameRl=[nm(1:end-2) stackSuffixU rlSuffix];

% change .mrc to .mrcs
if ~exist(sDatNameRl,'file')
    if exist(sDatName,'file')
        system(['mv ' sDatName ' ' sDatNameRl]);
     else
        error([sDatName ' not found']);
    end;
end
if unsubToo && ~exist(uDatNameRl,'file')
    if exist(uDatName,'file')
        system(['mv ' uDatName ' ' uDatNameRl]);
else
        error([uDatName ' not found']);
    end;
end

%%  
minSD=1.115;
maxSD=1.17;
if checkSD
    disp(['Reading ' sDatNameRl]);
    sImgs=ReadMRC(sDatNameRl);
%%
    % Find images having sufficient variance
    [n, n1, nim]=size(sImgs);
    sds=std(reshape(sImgs,n*n1, nim));
    subplot(2,1,1);
    plot(sds);
    ylabel('Image std');
    xlabel('Particle number');
    sdLower=sds(:)>minSD;
    disp([num2str(sum(~sdLower)) ' below minSD']);
    sdUpper=sds(:)<maxSD;
    disp([num2str(sum(~sdUpper)) ' above maxSD']);
    allActive=sdLower & sdUpper & activeFlags;
    disp([num2str(sum(allActive)) ' images active and in bounds.']);
%    clear('sImgs');  % Clear out the big array.
    
else
    allActive=activeFlags;
end;
return
%%
% Output file name
starName=[nm '.star'];
disp(['Output file: ' starName]);
%

fi=fopen(starName,'w');
% fi=1;

fprintf(fi,'%s %s\n#\n','# Original stack name: ',siName);
fprintf(fi,'data_%s\n',dataName);

fprintf(fi,'loop_\n');
fprintf(fi,'_rlnImageName\n');
fprintf(fi,'_rlnMicrographName\n');
fprintf(fi,'_rlnDefocusU\n');
fprintf(fi,'_rlnDefocusV\n');
fprintf(fi,'_rlnDefocusAngle\n');
fprintf(fi,'_rlnVoltage\n');
fprintf(fi,'_rlnAmplitudeContrast\n');
fprintf(fi,'_rlnSphericalAberration\n');
fprintf(fi,'_rlnMagnification\n');
fprintf(fi,'_rlnDetectorPixelSize\n');
fprintf(fi,'_rlnCtfFigureOfMerit\n');
fprintf(fi,'_rlnGroupNumber\n');
if unsubToo
   fprintf(fi,'_rlnReconstructImageName\n');
end;

%allActive=activeFlags & sdOk;
numActiveImages=sum(allActive);
defoci=NaN*ones(nim,1,'single');
for i=1:nim
    if allActive(i) % some signal
        imgName=sprintf('%05d@%s',i,sDatNameRl);
        imgUName=sprintf('%05d@%s',i,uDatNameRl);
        mi=si.mi{si.miIndex(i)};
        if any(setWeights)
            mi.weights=setWeights;
        end;
        ds=si.pixA/mi.pixA;
        detPixel=detPixel0*ds;  % we scale up the pixel size for downsampling.
        mag=detPixel*1e4/si.pixA;
        mName=[mi.baseFilename 'm.mrc'];
        ctf=mi.ctf(1);
        defocus=ctf.defocus;
        defoci(i)=defocus;
        deltadef=ctf.deltadef;        
        defU=round((defocus+deltadef)*1e4);
        defV=round((defocus-deltadef)*1e4);
        ang=mod(ctf.theta*180/pi,180);
        kV=mi.kV;
        if ampContrast>0  % force a new value
            alpha=ampContrast;
        else
            alpha=ctf.alpha;
        end;
        Cs=ctf.Cs;
        fom=si.sVesicle(i)*100;  % figure of merit is vesicle amplitude *100
        gn=si.miIndex(i);        % group number
        if unsubToo
            fprintf(fi,'%s %s %d %d %6.2f %g %5.3f %g %g %g %g %g %s\n',imgName,mName,defU,defV,ang,kV,alpha,Cs,mag,detPixel,fom,gn,imgUName);
        else
            fprintf(fi,'%s %s %d %d %6.2f %g %5.3f %g %g %g %g %g\n',imgName,mName,defU,defV,ang,kV,alpha,Cs,mag,detPixel,fom,gn);
        end;
    end;
end;
fprintf(fi,'\n');
fclose(fi);
subplot(2,1,2);
plot(defoci);
ylabel('defocus');
xlabel('Particle number');
disp(' written.');