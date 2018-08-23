% ConvertMCToMicrographs.m
% Assuming we're in the experiment directory, with allNames containing the
% mi file names (with path), read each dose-weighted .mrc file
% written by MotionCor2, pad
% and rotate it to conform to our Micrograph file standard.  Put them into
% a Micrograph_MC/ folder.  Note that MotionCor2's dose weighting reduces
% the amplitude and produces a HF rolloff of the power spectrum--we do not
% compensate for this.

numImgs=2;  % convert both images
overwrite=1;

nim=numel(allNames);
outDir='Micrograph_MC/';
CheckAndMakeDir(outDir);

for i=1:nim
    name=allNames{i};
    if exist(name,'file')
        disp(['mi file: ' name]);
        mi=ReadMiFile(name);
        for j=1:min(numImgs,numel(mi.imageFilenames))
            
            [pa,nm,ex]=fileparts(mi.movieFilename{j});
            SumName=[mi.moviePath 'Sum_' nm '_DW.mrc'];
%            SumName=[mi.moviePath 'Sum' nm '_DW.mrc'];
            if exist(SumName,'file')
                outName=[outDir mi.imageFilenames{j}];
                disp(['      ' outName]);
                if overwrite || ~exist(outName,'file');
                    
                    m0=ReadMRC(SumName);
                    m=rot90(Crop(m0,mi.imageSize,0,mean(m0(:))),3);
                    imags(BinImage(m,4));
                    title(SumName,'interpreter','none');
                    drawnow;
                    
                    WriteMRC(m,mi.pixA,outName);
                else
                    disp('      --file exists, not overwritten.');
                end;
            else
                disp(['      --input file not found: ' SumName]);
            end;
        end;
    end;
end;