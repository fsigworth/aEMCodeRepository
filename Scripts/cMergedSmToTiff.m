% scMergedSmToTiff.m
% Read some .mrc files in Merged_sm and write them out as .tif images.
% 
inDir='Merged_sm/';
outDir='Merged_tif/';
outInfoDir='Info_tif/';
inExt='.mrc';
outExt='.tif';

goodIDs= ...
   {'103_0001'
    '103_0002'
    '104_0001'
    '104_0002'
    '104_0003'
    '105_0001'
    '105_0002'
    '106_0001'};
basePtr=34; % hardwired cutoff of name


d=dir(inDir);
nim=0;
goodNames=cell(0,1);
goodBaseNames=cell(0,1);
for i=1:numel(d)
    q=contains(d(i).name,goodIDs);
    if any(q) % we found something
        name=d(i).name;
        [pa,nm,ext]=fileparts(name);
        if ~strcmp(ext,inExt) 
            disp([name ' skipped.']);
            continue;
        end;
        inName=[inDir name];
        nim=nim+1;
        [~,goodBaseNames{nim,1}]=fileparts(inName(1:basePtr));
%         goodBaseNames{nim,1}=inName(1:basePtr);
        goodNames{nim,1}=inName;
        outName=[outDir nm outExt];
        m=ReadMRC(inName);
        WriteJpeg(m,outName);
%         disp([inName ' > ' outName]);
    end;
end;
goodBaseNames
%%
vm=[0 0 .8 0 0];
ds=4;
for i=1:nim-1
    miName=['Info/' goodBaseNames{i} 'mi.txt'];
    mi=ReadMiFile(miName);
    mv=numel(mi.vesicle.x);
     mi.vesicle.s=ones(mv,1);
     mi.vesicleModel=vm;
    vImg=-meMakeModelVesicles(mi,mi.padImageSize/ds,0,0,0,1);
    imags(vImg);
    outName=[outDir goodBaseNames{i} 'v.tif'];
    title(outName,'interpreter','none');
    drawnow;
    WriteJpeg(vImg,outName);
    WriteMiFile(mi,[outInfoDir goodBaseNames{i} 'mi.txt']);
end;




