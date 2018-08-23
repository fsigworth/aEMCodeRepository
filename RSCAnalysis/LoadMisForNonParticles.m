% LoadMisForNonParticles
% Create an si file containing the mi structures from selected files

fInds=[1:10 51:60];
cd /Users/fred/EMWork/Hideki/151117/KvLipo80slot3
d=dir('Info');
d(1:4)=[];  % ignore the first ones.
si=struct;
si.mi=cell(0,1);
si.miIndex=zeros(0,1,'int16');
nut=0;
ind=1;
for i=1:numel(fInds)
    ind=fInds(i);
    name=['Info/' d(ind).name];
    mi=ReadMiFile(name);
    nu=sum(mi.particle.picks(:,3)==48);
    disp([name '  ' num2str(nu) ' blanks']);
    si.mi{i,1}=mi;
    si.miIndex(ind:ind+nu-1)=1;
    ind=ind+nu;
end;
disp([num2str(ind-1) ' total ''particles''.']);
save Stack/NonPartsi.mat si
