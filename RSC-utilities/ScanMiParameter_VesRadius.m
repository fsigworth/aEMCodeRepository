% ScanMiParameter_VesRadius

% field='frameDose';
% field='Defocus, µm';
field='vesicle.r';
nv=1;  % pick up this many values
xbins=0:10:990;  % histogram bins
nbins=numel(xbins);

% Have the user select some mi files
[fname, pa]=uigetfile('*mi.txt','Select mi files to examine','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
%%
cd(pa);
nim=numel(fname);
val=zeros(nim,nbins);
figure(1);
subplot(111);
mis=cell(0,1);
for i=1:nim
    %     load(fname{i});
    mi=ReadMiFile(fname{i});
    mis{i,1}=mi;
    if numel(mi.vesicle.x)>0
            okVes=all(mi.vesicle.ok,2);
            radii=mi.vesicle.r(:,1)*mi.pixA;
%             radii=mi.vesicle.r(:,1);
            val(i,:)=hist(radii,xbins);
            
    bar(xbins,val(i,:));
    title(fname{i},'interpreter','none');
    drawnow;
    end;
end;

bar(xbins,sum(val,1));
xlabel('Vesicle radius, Å');
ylabel('Frequency');

save AllMis.mat mis
