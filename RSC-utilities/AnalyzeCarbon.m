[fname, pa]=uigetfile('*.*','Select image files','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;

cd(pa);

nFiles=numel(fname);
for i=1:nFiles
    [m pixA]=ReadEMFile(fname{i});
    ms(:,:,i)=RemoveOutliers(m);
    imacs(ms(:,:,i));
    title(fname{i},'interpreter','none');
    drawnow;
end;

%% Prewhiten
mw=mePreWhiten(ms);

%%

q=mw(3001:4024,501:1524,1);
imacs(q)
q1=Crop(q,768);
imacs(q1)
%%
mg=RemoveGoldParticles(mw);




%%
semilogy(RadialPowerSpectrum(q1));
