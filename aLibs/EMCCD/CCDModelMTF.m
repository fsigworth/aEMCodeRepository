function mtf=CCDModelMTF(f,iCamera)
% function mtf=CCDModelMTF(f)
% Return the analytical fit to the MTF of our CCD at 200kV
% as determined with MeasureMTF3.
% f is a frequency vector, with 0.5 (maximum value) being Nyquist.
% iCamera=1 : Yale F20 US4000
% iCamera=2 : Okazaki JEM 2200 TVIPS camera
% iCamera=3 : DE12
%
% f is a frequency vector of any dimension, with 0.5 (maximum value)
% being Nyquist.

if nargin<2
    iCamera=1;
end;

switch iCamera
    case 1
        f1=.0857;
        f2=.2250;
        f3=.3874;
        f=abs(f);
        mtf=(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));
        
    case 3
        f0=0.01513;
        a0=0.45748;
        f1=0.20599;
        f2=0.38772;
        f3=0.51004;
        mtf=a0./(1+(f/f0).^2)+(1-a0)*(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));

    case 4
        p=[0.026528     -0.21526     0.043531     0.099322      0.27167];
        f0=p(1); a0=p(2); f1=p(3); f2=p(4); f3=p(5);
        mtf=a0./(1+(f/f0).^2)+(1-a0)*(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));
    
    case 5 % k2 camera, fit to 7e/pix/s curve in Ruskin et al. JSB13
        x=f*2;  % fraction of Nyquist.  Fit from Curve Fitting Tool.
%         The value of mtf at f=0 is 0.772; includes counting loss.
       p1 =      0.5704; %  (0.4634, 0.6773)
       p2 =      -1.312; %  (-1.454, -1.17)
       p3 =      0.9687; %  (0.8886, 1.049)
       q1 =       -1.64; %  (-1.783, -1.498)
       q2 =       1.255; %  (1.152, 1.359)
        mtf = (p1*x.^2 + p2*x + p3) ./ (x.^2 + q1*x + q2);
        
    case 6
        p = [0.021788     0.021042     0.063745      0.11211      0.24927];
    f0=p(1); a0=p(2); f1=p(3); f2=p(4); f3=p(5);
    mtf=a0./(1+(f/f0).^2)+(1-a0)*(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));
    otherwise
        error(['Invalid camera index ' num2str(iCamera)]);
end;

