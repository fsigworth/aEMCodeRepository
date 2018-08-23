% TestVesicleDistributions
%% 
[fnames pa]=uigetfile('*mi.mat','Select an mi file','multiselect','on');
fnameCell=cell(1);
if isa(fnames,'cell')
    nf=numel(fnames);
    fnameCell=fnames;
elseif isnumeric(fnames)
    return
else
    nf=1;
    fnameCell{1}=fnames;
end;
%% 

xs=10:10:400;
nx=numel(xs);
h=zeros(nx,1);
    cd(pa);
    for i=1:nf
    load([AddSlash(pa) fnameCell{i}]);
    disp(fnameCell{i});
    radii=mi.pixA*mi.vesicle.r;
    h=h+hist(radii,xs)';
    end;
bar(2*xs,h);
title(pa,'interpreter','none');
xlabel('Vesicle diameter, Å');

return

%% Create the figure

% TestVesicleDistributions
h35=h;
names35=fnameCell;
% TestVesicleDistributions


h20=h;
names20=fnameCell;
sum(h)
sum(h35)
h351=975/458*h35;

%%
bar(2*xs,[h20 h351]);
xlabel('vesicle diameter, Å');
ylabel('frequency');
legend('20PE','35PE');
text(400,140,'PE-20');
text(400,120,names20,'interpreter','none');
text(560,140,'PE-35');
text(560,120,names35,'interpreter','none');