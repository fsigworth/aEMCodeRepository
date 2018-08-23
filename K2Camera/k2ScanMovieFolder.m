% k2ScanMovieFolder

% TestAnalyzeK2Movies

cd('/Volumes/TetraSSD/130822/GluR2GFP_glutamate_slot4/movie_frames/');
d=dir;
xs=0:200;
for i=1:numel(dir)
    [pa nm ex]=fileparts(d(i).name);
    if strcmp(ex,'.mrc')
        disp(d(i).name);
        m=ReadMRC(d(i).name);
        h=hist(mf(:),xs);
        figure(1);
        bar(xs,h);
        drawnow;
        hs=sum(h);
        %%
        cutoff=find(h<hs*1e-7,1,'first')
        mf(mf>=cutoff)=0;
        mfs=sum(mf,3);
        figure(2);
        SetGrayscale;
        imacs(mfs);
        drawnow;
        disp(' ');
    end;
end;

%