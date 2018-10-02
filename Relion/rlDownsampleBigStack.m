% rlDownsampleBigStack.m
scale=-1;
ds=4;
d=dir;

for ind=13:14
    stackName=d(ind).name;
    disp(stackName);
    [m,s]=ReadMRC(stackName,1,1);
    n=size(m);
    
    nim=s.nz;
    name=d(ind).name
    block=5000;
    nBlocks=ceil(nim/block)
    imgs=[];
    
    for i=1:nBlocks
        iStart=(i-1)*block+1;
        disp(iStart);
        iNum=min(block,nim-iStart+1);
        m0=ReadMRC(d(ind).name,iStart,iNum);
        m1=scale*Downsample(m0,n/ds,1);
        if i==1
            imgs=m1;
        else
            imgs(:,:,end+1:end+iNum)=m1;
        end;
    end;
    flag='';
    if ind==14
        flag='u';
    end;
    outName=['ds128t' flag 'stack.mrc'];
    disp(outName);
    WriteMRC(imgs,s.pixA*ds,outName);
end;
%%
load(d(12).name);  % load the si.mat file
%
si=rsStackDownsample(si,[],128);
save('dssi.mat','si');
disp('dssi.mat');
