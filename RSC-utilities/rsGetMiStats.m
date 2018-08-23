% rsGetMiStats.m
% from a text file of mi file names, extract values of a field (in this
% case, mi.mergeMatrix,



[name, pa]=uigetfile('*.*','Select a text file');
    if isnumeric(pa) % File selection was cancelled
        return
    end;
f=fopen([AddSlash(pa) name]);
data=zeros(3,3,0);
k=0;
if f>0
    while ~feof(f)
        line=fgetl(f);
        if exist(line,'file')
            load(line);
            if numel(mi.mergeMatrix)>9
                for j=2:size(mi.mergeMatrix,3)
                    k=k+1;
                    data(:,:,k)=mi.mergeMatrix(:,:,j);
                    if mod(k,1000)==0
                        disp(k);
                    end;
                end;
            end;
        end;
    end;
else
    error('text file couldn''t be opened');
end
%%
subplot(2,2,1);
bins=.91:.002:1.19;
h=hist(squeeze(data(1,1,:)),bins);
bar(bins,h);
title('T(1,1)')

subplot(2,2,3);
plot(bins,cumsum(h));

subplot(2,2,2);
bins=-.03:.001:.03;
h=hist(squeeze(data(1,2,:)),bins);
bar(bins,h);
title('T(1,2)');

subplot(2,2,4);
plot(bins,cumsum(h));