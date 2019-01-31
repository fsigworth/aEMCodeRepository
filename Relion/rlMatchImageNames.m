% rlMatchImageNames(st1,st0)
st1=d1.rlnImageName;
n1=numel(st1);
st0=d.rlnImageName;

[sst0,sp0]=sort(st0);
[sst1,sp1]=sort(st1);
p1=1; % found pointer
fp1=zeros(n1,1);
for i=1:4
    p1=find(sst0(p1:end),sst1{i},'first');
    if numel(p1)>0
        fp1(i)=p1;
        p1
    end;
end;

%%
st1=string(d1.rlnImageName);
n1=numel(st1);
st0=string(d.rlnImageName);
pt=zeros(n1,1);
for i=1:n1
    pt(i)=find(strcmp(st1(i),st0),1);
end;


