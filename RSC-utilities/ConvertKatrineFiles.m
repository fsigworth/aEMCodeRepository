% ConvertKatrineFiles.m
% Convert Katrine's .mat files to mrc vesicle model images.
doWrite=1;

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
dir1='Vesicle_models2/';  % Katrine's data
dir2='Vesicles2/';        % Where we'll write .mrc images and .mat logicals
dir1='Vesicle_models_R8A3PC24/';
dir2='VesiclesR8A3/';
CheckAndMakeDir(dir2,1);
dirj2=[dir2 'jpeg/'];
CheckAndMakeDir(dirj2,1);

d1=dir(dir1);

n1=numel(d1);
for i=1:n1
    kName=d1(i).name;
    [pa,nm,ex]=fileparts(kName);
    if strcmp(ex,'.mat')  % it's a mat file
        s=load([dir1 kName]);
        mv=s.mMicrographVesiclesDifference';
        ms=s.mMicrographVesiclesSubtracted';
        outName=[dir2 nm '.mrc'];
%         imaga(256*BinImage(mv,4)+150);
%         title(outName,'interpreter','none');
%         drawnow;
        disp(outName);
        if doWrite
%             WriteMRC(mv,0,outName);
        outJName0=[dirj2 nm '.jpg'];
        outJName1=[dirj2 nm 'm.jpg'];
        outJName2=[dirj2 nm 's.jpg'];
            WriteJpeg(BinImage(mv,2),outJName0);
            imaga(150*BinImage(ms,2)+150);
            drawnow;
%             WriteJpeg(BinImage(ms+mv,2),outJName1);
%             WriteJpeg(BinImage(ms,2),outJName2);
        end;
    end;
end;
 