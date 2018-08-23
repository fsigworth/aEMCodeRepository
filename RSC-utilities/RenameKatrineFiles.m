% RenameKatrineFiles
% Match up Katrine's vesicle models with the original merged images, and
% write out vesicle model images with canonical filenames.

pattern='070215';  % End of fixed part of Katrine's filenames

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
dir0='Merged/';          % merged images to match
dir1='Vesicle_models/';  % Katrine's data
dir2='Vesicles/';        % Where we'll write .mrc images and .mat logicals

d0=dir(dir0);
d1=dir(dir1);

n1=numel(d1);
n0=numel(d0);
serno=zeros(n1,1);
kNames=cell(n1,1);
for i=1:n1
    kName=d1(i).name;
    kNames{i}=kName;
    p=strfind(kName,pattern);
    q=strfind(kName,'.mat');
    if numel(p)>0 && numel(q)>0
        p=p(1)+numel(pattern);
        q=q(1)-1;
        serno(i)=str2double(kName(p:q));
    end;
end;
[vals,inds]=sort(serno);
kNames=kNames(inds);
zers=(vals==0);
inds(zers)=[];
% vals(zers)=[];
kNames(zers)=[];
%%
vImgs=zeros(480,480,'single');
mImgs=zeros(480,480,'single');

for j=1:max(n0,n1)
    if j<=numel(kNames)
        s=load([dir1 kNames{j}]);
        vImg=BinImage(s.mMicrographVesiclesDifference,4)';
        vImg=vImg-mean(vImg(:));
        vImgs(:,:,j)=vImg;
    end;
    if j<=n0-3
        mImg=BinImage(ReadZTiff([dir0 d0(j+3).name]),4);
        mImg=mImg-mean(mImg(:));
        mImgs(:,:,j)=mImg;
    end;
    %     subplot(221);
    %     imags(mImg);
    %     title(d1(i).name);
    %     subplot(222);
    %     imags(vImg);
    %     title(d0(j+3).name);
    %     pause;
    cc(j)=(vImg(:)'*mImg(:))/(sqrt(vImg(:)'*vImg(:))*sqrt(mImg(:)'*mImg(:)));
    disp([j cc(j)]);
end;
%%
np=480^2;
nim0=size(mImgs,3);
nim1=size(vImgs,3);
vIVec=reshape(vImgs,np,nim1);
mIVec=reshape(mImgs,np,nim0);
ccM=vIVec'*mIVec;  % rows are vImg, cols are mImg
imagesc(ccM);
axis ij
ylabel('Ves image number');
xlabel('Merged image number');

%%
[vals,matches]=max(ccM);
matches(vals<1000)=0;  % matches are the best vImg for each mImg

for j=1:n0-3
    mName=d0(j+3).name;
    [pa,nm,ex]=fileparts(mName);
    nm(end)=[];
    nm(end)='v';  % replace the last character
    vName=[nm '.mrc'];
    okName=[nm(1:end-1) 'ok.mat'];
    if matches(j)>0
        oldVName=kNames{matches(j)};
        disp([mName '   ' oldVName '   ' vName]);
        s=load([dir1 oldVName]);
        vImg=s.mMicrographVesiclesDifference';
        vOk=~logical(s.vIdxDiscarded);
%         subplot(221);
%         imags(vImg);
%         subplot(222);
%         [mImg,s]=ReadZTiff([dir0 mName]);
%         mImg=mImg-mean(mImg(:));
%         imags(mImg);
%         subplot(223);
%         imags(mImg-vImg);
      drawnow;
        save([dir2 okName],'vOk');
%         WriteMRC(vImg,1,[dir2 vName]);

    end;
    
end;
