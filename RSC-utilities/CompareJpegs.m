% CompareJpegs
% Utility to manually compare vesicle subtraction making use of jpeg images
% in a set of directories.  For example I ran FindVesicle with four
% different minimum vesicle amplitudes (findVesicleAmps) using
% four different working directories (findVesicleDirs).  We display 



defaultStartIndex=250;
basePath='/ysm-gpfs/scratch60/fjs2/160909/';
findVesicleAmps=[5e-4 6e-4 7e-4 8e-4];
findVesicleDirs={'KvLipo121_2w10_v3.5/';'KvLipo121_2w10_v3.6/';
    'KvLipo121_2w10_v3.7/';'KvLipo121_2w10_v3.8/'};

cd([basePath findVesicleDirs{1}]);
d=dir('Jpeg');
names=cell(0,1);
for i=3:numel(d)
    names{i-2}=d(i).name;
end;
nf=numel(names);
nj=numel(findVesicleDirs);
indFileName=[basePath 'inds.mat'];
if ~exist(indFileName,'file')
    disp(['Creating the index array ' indFileName]);
inds=ones(nf,1);
else
    load(indFileName);
inds=inds(:,1);
    save([basePath 'indsBak.mat'],'inds'); % save a backup copy
end;
startIndex=find(inds==0,1)
if numel(startIndex)==0
    startIndex=defaultStartIndex;
end;
startIndex=startIndex(1);

defaultIndex=1;

figure(2);
i=startIndex;
done=false;
while ~done
    for j=1:nj
        cd([basePath findVesicleDirs{j} 'Jpeg/']);
        if exist(names{i},'file')
            img=imread(names{i});
        else
        img=0;
        end;
        mysubplot(2,2,j)
        image(img);
        title([findVesicleDirs{j} names{i}],'interpreter','none');
        axis off;
    end;
    
    q=MyInput('Index ', defaultIndex);
    w=q(1);
    if q<0
        i=max(1,i-1);
    elseif q==0
        inds(i)=0;
        break
    else
        inds(i)=q;
        defaultIndex=q;
        i=i+1;
    end
    if i>nf
        done=false;
    end;    
end;
%%
disp(['Saving ' indFileName]);
save(indFileName,'inds');
