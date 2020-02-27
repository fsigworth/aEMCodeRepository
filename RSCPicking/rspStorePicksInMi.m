function mi=rspStorePicksInMi(mi,picks,ptrs,doDisplay)
if nargin<4
    doDisplay=1;
end;
% 
% particle types are the following
% 0-1 not displayed
%  0 = null
%  1 = deleted from auto
% 2 = marked vesicle
% 3 = bad vesicle
% 16 manual picking
% 32 auto picking
% 48 background
% 
% We  store all non-null entries in mi.particle.picks
% We also mark mi.vesicle.ok(:,4) according to bad vesicles.

q=picks(:,:,3)>0;  % Find all the nontrivial entries, with type ~=0.
nk=sum(q(:));      % number of nontrivial entries
ne=size(picks,3);  % Number of elements in each row of picks (9?)
mi.particle.picks=single(zeros(nk,ne));
mi.particle.ok(:,4)=true;  % default is all are ok.
k=0;
for j=1:numel(ptrs)  % There are 7 elements of ptrs, 7 classes of picks
    for i=1:ptrs(j)
        c=picks(j,i,:);
        c=c(:)';  % make into a row vector
        if c(3)>0 % nontrivial entry
            k=k+1;
            mi.particle.picks(k,:)=c;
            if c(3)==3 % bad vesicle
                mi.vesicle.ok(c(4),4)=false;
%                 --this actually does nothing as ok(:,4) is usually 0
%                 anyway.
            end;
        end;
    end;
end;
if doDisplay
    [np,npv,npa]=rspGetParticleStats(mi);
    sp='•  ';
      disp([sp num2str(np) ' particles']);
        disp([sp num2str(npv,3) ' particles/vesicle']);
        disp([sp num2str(npa,3), ' particles/area']);
end;
