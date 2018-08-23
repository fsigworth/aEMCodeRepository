% k2MoviesToJpegs
% Simple script to read movies and deposit their sums as Jpegs

fixedGrayscale=1;
add=-160;
scl=6;


% Have the user select some movie files
[names, pathName]=uigetfile('*.*','Select movie files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)
    names={names};
end;
cd(pathName); % make a jpeg dir inside the current dir
jpegDir='Jpeg/';
ok=CheckAndMakeDir(jpegDir,1);
figure(1);
clf;
nim=numel(names);
%%
for i=1:nim
    name=names{i};
    disp(['Reading ' name]);
    [m,s,ok]=ReadMovie(name);
    if ok
        [pa,nm,ex]=fileparts(name);
        jname=[jpegDir nm '.jpg'];
        disp(['Writing ' jname]);
        img=BinImage(sum(single(m),3),2);
        imgMean=mean(img(:))
        %%
        if fixedGrayscale
            ms=uint8(rot90(img*scl+add));
            imwrite(ms,jname);
            imaga(ms);
        else % use our autoscaling
            WriteJpeg(img,jname);
            imaga(imscale(BinImage(img,2),256,1e-3));
        end;
        title(jname,'interpreter','none');
        axis off;
        drawnow;
    end;
    
end;
disp('Done.');

