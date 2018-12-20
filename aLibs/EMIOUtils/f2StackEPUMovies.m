% f2StackEPUMovies

extension='.mrc';
% [nm,dataPath]=uigetfile(extension,'Select the first movie file');
 d=dir;
outPath='/Users/fred/movies/';
CheckAndMakeDir(outPath,1);

fi=1;
while fi<numel(d)
    [pa,nm,ex]=fileparts(d(fi).name);
    
    if ~d(fi).isdir && strcmp(ex,extension) && strndcmp(nm,'n0',2)
        %         We've found the first of a set of images
        baseName=nm(1:end-10);
        disp(baseName);
        [m,s]=ReadMRC(d(fi).name);
        fprintf('.');
        ni=1; % image counter
        nj=1; % file counter
        nextName=d(fi+nj).name;
        while fi+nj < numel(d) && strncmp(nextName,baseName,numel(baseName))
            [pa,nm,ex]=fileparts(nextName);
            if strcmp(ex,extension)
                m(:,:,ni+1)=ReadEMFile(d(fi+nj).name);
                fprintf('.');
                ni=ni+1;
            end;
            nj=nj+1;
            nextName=d(fi+nj).name;
        end;
        fprintf('\n');
        ms=sum(single(m),3);
        msm=Downsample((ms),size(ms)/4);
        imags(msm);
        title(baseName);
        drawnow;
%         WriteMRC(m,0,[outPath baseName '.mrcs'],1);
        WriteMRC(msm,0,[outPath baseName 'sm.mrc']);
        WriteJpeg(msm,[outPath baseName 'sm.jpg'],1e-5);

        disp([num2str(ni) ' frames.']);
        fi=fi+ni;
    else
        fi=fi+1;
    end;
end;

