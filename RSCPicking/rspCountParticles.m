function counts=rspCountParticles(picks)
% returns the vector of the number of (manually, auto, vesicles, bkd) picked
% particles.
% 
% particle types are the following
% 0-1 not displayed
%  0 = null
%  1 = deleted from auto
%  2 = vesicles
%  3 = bad vesicles
% 16-31 manual picking (we use only 16 and 17 at present)
% 32-47 auto picking (we use only 32 at present)
% 48-63 background (we use only 48 at present)
flags=[16 17 32 2 48];
temp=0*flags;
for i=1:numel(flags)
    q=picks(:,:,3)==flags(i);
    temp(i)=sum(q(:));
end;
counts=[temp(1)+temp(2) temp(3) temp(4) temp(5)];

