% rsAssignGroups.m

% Assign groups according to intensity measures and defocus values.

% cd Picking_9
% load allMis_holemasked_i2.mat
% load defJumpStarts.mat

nim=numel(allMis);
idat=zeros(nim,3);
for i=1:nim
    mi=allMis{i};
    idat(i,:)=mi.ok([9 10 14]); % fraction, ves intensity, top 1% intensity
end;

figure(1);
subplot(311);
hist(idat(:,1),200); % fractions
subplot(312)
hist(idat(:,2:3),200); % intensities

sel=true(nim,2);
d2=.18; % approx displacement between vesicle and max intensities
sel(:,1)=idat(:,1)>.1; % vesicle area is >10%
sel(:,2)=idat(:,2)>7.5-d2 & idat(:,2)<8.5-d2;  % vesicle intensity is reasonable
sel(:,3)=idat(:,3)>7.5 & idat(:,3)<8.5;  % max intensity is reasonable

id2=idat(all(sel(:,1:2),2),2);
id3=idat(sel(:,3),3);

subplot(312)
hist(id2,200);
subplot(313)
hist(id3,200);

id2Shift=median(id3)-median(id2); % precise value of d2
id2=id2+id2Shift;

[h2,bins]=hist(id2,200);
h3=hist(id3,bins);

bar(bins,[h2' h3']);

% compute the composite intensity as the mean of the peak and vesicle
% intensities, when vesicle intensities are ok.

peakIntensOk=sel(:,3);
vesIntensOk=all(sel,2);

compIntens=mean([idat(:,2)+id2Shift idat(:,3)],2); % assign the mean values
compIntens(~vesIntensOk)=idat(~vesIntensOk,3); % pick up the peak value when vesicle data are missing.
compIntens(~peakIntensOk)=NaN;

hist([compIntens(peakIntensOk) idat(peakIntensOk,3)],200);

for i=1:nim
    allMis{i}.ok(8)=compIntens(i); % This is where we put the composite intensity.
    allMis{i}.ok(7)=any(i==defJumpStarts); % =1 if a defocus jump starts after this.
end;

% save allMis_holemasked_i2.mat



%% Assign class numbers

nmi=numel(allMis);
vals=zeros(nmi,15);
for i=1:nmi
    z=allMis{i}.ok(1:15);
    vals(i,:)=z;
end;

% valid intensity, resolution and defocus.
allBads=isnan(vals(:,8)) & (res>7) & defs>3.6;

%  Set the selection criteria here:
q=1./(1+vals(:,1).^2); % switch function favoring overlap <=1
% vals(:,8) is the composite intensity
selVal=vals(:,8)+(8.4-vals(:,8)).*q/2;
selVal2=vals(:,1);
figure(2);
clf;
hist(selVal,200);

figure(3);
selVal2(selVal2>20)=NaN;
plot(selVal,selVal2,'.');

inds=find(selVal>8.2 & selVal<8.5 & defs<3.6 & ~isnan(selVal) & res<7);

% we'll sort into four classes for defocus and intensity.
defBounds=[-1 1 1.8 2.3 2.7 3.6];
intBounds=[-1 7.6 7.9 8.1 8.2 8.5];

defClass=zeros(nim,1);
intClass=zeros(nim,1);
df=defs;
df(isnan(df))=0;
sv=selVal;
sv(isnan(sv))=0;
for i=1:nim
    defClass(i)=find(df(i)>defBounds,1,'last')-1;
    intClass(i)=find(sv(i)>intBounds,1,'last')-1;
end;
    defClass(defClass>4)=0;
    intClass(intClass>4)=0;

classNumber=(defClass>0).*(intClass>0).*(defClass*4+intClass-4);
classMatrix=zeros(nim,4,4);
for i=1:nim
    if intClass(i)>0 & defClass(i)>0
        classMatrix(i,intClass(i),defClass(i))=1;
    end;
end;

for i=1:nim
    allMis{i}.ok(20)=classNumber(i);
end;

save Picking_9/allMis_holes_i2_ov_cls.mat
