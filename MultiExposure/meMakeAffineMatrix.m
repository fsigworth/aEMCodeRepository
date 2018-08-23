function T=meMakeAffineMatrix( P )
% Construct a 3x3 affine transform matrix from the parameters in the
% vector P, with elements
% 1 rotation angle (radians)
% 2 magnification
% 3 x-stretch
% 4 xy-stretch
% 5 x shift
% 6 yshift
% Defaults for all parameters are 0 except magnification=1.
ne=numel(P);
if ne<2
    P(2)=1;
    ne=2;
end;
if ne<6
    P(ne+1:6)=0;
end;

c=cos(P(1));
s=sin(P(1));
T=zeros(3,3);
mag=P(2);
sxx=P(3);
sxy=P(4);
T=P(2)*[c -s  0 ; s c  0 ; 0 0 0]...
    +[sxx sxy P(5); sxy -sxx P(6); 0 0 1];
