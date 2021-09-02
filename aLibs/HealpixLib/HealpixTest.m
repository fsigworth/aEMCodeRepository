% Healpix Test
% Healpix order, number of points:
%      1    12
%      2    48 = 4*12
%      3   108
%      4   192 = 4*48
%      5   300
%      6   432
%      7   588
%      8   768 = 4*192
%      9   972
%     10  1200
%     11  1452
%     12  1728
%     13  2028
%     14  2352
%     15  2700
%     16  3072 = 4*768
        

n = 4;
% n = 4;
% n = 16;
for n=1:16
I = HealpixGenerateSampling(n, 'rindex');
S = HealpixGenerateSampling(n, 'scoord');

% S= [theta phi] in radians.

C = SphToCart(S);
subplot(121);
plot3(C(:,1),C(:,2),C(:,3), '.')
%plot3(C(:,1),C(:,2),C(:,3))
axis equal

phis=S(:,2)*180/pi;
thetas=S(:,1)*180/pi;
subplot(122);

plot(phis,thetas,'o');

title(size(S,1));
drawnow;

np=size(S,1);
disp([n np]);
end;