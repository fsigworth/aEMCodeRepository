% rlTrimBadNorms

tolerance=0.2;
% inDir='/Select/job041/';
% cd(inDir);
inName='particles.star';
inFile=inName;
% outDir=inDir;
outName='particles_trimmed2.star';
outFile=outName;

% [nms,dats]=ReadStarFile(inFile);
d=dats{2};
med=median(d.rlnNormCorrection);
active=abs(d.rlnNormCorrection/med-1)<tolerance;
disp(['Rejected ' num2str(sum(~active)) ' of ' num2str(numel(active))]);

ok=MyInput('Proceed? ','y');
if ok=='y'
    if isa(dats{1}.rlnOpticsGroupName,'char')
        dats{1}.rlnOpticsGroupName={dats{1}.rlnOpticsGroupName};
    end;
    disp(['Writing ' outFile]);
    WriteStarFileStruct(dats,nms,outFile,{[] active});
    disp('done.');
else
    disp('Nothing written.');
end;
