% Analyze vesicle models

% [mi,name]=ReadMiFile;
%% Alternative
cd('/Users/fred/EMWork/Hideki/170810-2/Info');
name='vp_0001_Aug10_14.58.56mi.txt';
mi=ReadMiFile(name);

ind=22;
nx=196;

% squeeze(mi.vesicle.s(ind,:,:))
%%
v0=meMakeModelVesicles(mi,960,ind);
imags(v0)

x0=ExtractImage(v0,round([mi.vesicle.x(ind) mi.vesicle.y(ind)]/4),nx);
%%
for i=1:3
    mi1=mi;
    mi1.vesicle.s=mi1.vesicle.s*0;
    mi1.vesicle.s(:,:,i)=mi.vesicle.s(:,:,i);
    v1=meMakeModelVesicles(mi1,960,ind);
    imags(v1); drawnow;
    x(:,:,i)=ExtractImage(v1,round([mi.vesicle.x(ind) mi.vesicle.y(ind)]/4),nx);
    pause
end;
plot(sect(x));
figure(2);
imags(sum(x,3));
% 
% mi2=mi;
% mi2.vesicle.s(:,:,1)=0;
% v2=meMakeModelVesicles(mi2,960,22);
% imags(v2)

