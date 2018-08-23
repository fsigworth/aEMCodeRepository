% Analyze vesicle models 2
% try various exPos values

% [mi,name]=ReadMiFile;
%% Alternative
cd('/Users/fred/EMWork/Hideki/170810-2');

origName='Info/vp_0001_Aug10_14.58.56mi.txt';
miOrig=ReadMiFile(origName);
vFields=fieldnames(miOrig.vesicle);
mi0=miOrig;

% Select just some of the vesicles to fit.
inds=13:22;
for i=1:numel(vFields)
    if size(mi0.vesicle.(vFields{i}),1)>1
    mi0.vesicle.(vFields{i})=mi0.vesicle.(vFields{i})(inds,:,:);
    end;
end;

mi0.vesicleModel=[];

% name='Info/up_0001_ves22ami.txt';
% mi0=ReadMiFile(name);

positions=[30:5:45];
np=numel(positions);
vesFits=cell(0);
% errs=zeros(np,44);
% amps=zeros(np,44,2);

for i=1:np
    xName=['Info/up_0001_ves22a' num2str(i,'%02d') 'mi.txt']

    mi=mi0;
    WriteMiFile(mi,xName); % store the startin mi file by the new name.
%     mi.vesicle.extraPeaks=[-positions(i) positions(i)];
    pars=struct;
%     pars.peakPositionA=[-1 1]*positions(i); % centers of extra peaks, in angstrom
    pars.peakPositionA=[-37 0 positions(i)]; % centers of extra peaks, in angstrom
%     pars.rTerms=inf;
    pars
    rsRefineVesicleFits({xName},pars);  % This updates the mi file.
    mi=ReadMiFile(xName);
    disp(mi.vesicle.extraPeaks);
%     disp(squeeze(real(mi.vesicle.s(22,1,2:3)))');
    disp(squeeze(real(mi.vesicle.s(1,1,2:3)))');
    if i==1
        errs=mi.vesicle.err';
        amps=squeeze(abs(mi.vesicle.s(:,1,:)));
        vesFits={mi.vesicle};
    else
        errs(i,:)=mi.vesicle.err';
        amps(:,:,i)=squeeze(abs(mi.vesicle.s(:,1,:)));
        vesFits{i}=mi.vesicle;
    end;
    disp(['model: ' num2str(positions(i)) '  ' num2str(errs(i))]);
end;

figure(3);
plot(positions,sum(errs,2));
% squeeze(mi.vesicle.s(ind,:,:))
%%
