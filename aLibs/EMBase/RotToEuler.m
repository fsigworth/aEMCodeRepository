function ea=RotToEuler(m);
% Given the 3x3xni rotation matrices m, compute the 3xni Euler angles
% [phi theta psi]. Note that the result is non-unique when theta is 0 or
% 180, in which case we assign only psi to be nonzero.
% Based on the equations of Appendix A of Heymann et al. JSB 2005
thetaTol=.9999;
ni=size(m,3);
ea=zeros(3,ni);
singular=abs(m(3,3,:))>thetaTol;
ea(1,:)=atan2d(m(3,2,:),m(3,1,:));
ea(2,:)=acosd(m(3,3,:));
ea(3,:)=atan2d(m(2,3,:),-m(1,3,:));
% Handle singular case
ea(1:2,singular)=0;
ea(3,singular)=m(3,3,singular).*atan2d(-m(2,1,singular),m(1,1,singular));
% ea(3,singular)=atan2d(-m(2,1,singular),m(1,1,singular));

