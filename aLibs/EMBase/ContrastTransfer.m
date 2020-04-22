function [c, chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha)
% function [c, chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha)
% function [c, chi]=ContrastTransfer(s, CPars, doComplex)
% Compute the contrast transfer function corresponding
% to the given spatial frequency s (in A^-1); lambda is in A, Cs is in mm,
% B in A^2 and alpha in radians.  Defocus is in microns.
% The returned chi variable is -1/pi times the argument to the sin function.
% Hence CTF=sin(-pi*chi).*exp(-B*s.^2). this form is convenient because
% chi=-1 at the first zero, -2 at the second, etc. as long as Cs effects
% are negligible (otherwise it's non-monotonic).
% In the alternate form, CPars is a structure containing those parameters.
% CPars.lambda
% CPars.defocus
% CPars.Cs
% CPars.B
% CPars.alpha
%   and, optionally, to handle astigmatism,
% CPars.deltadef
% CPars.theta
%   and, optionally to handle phase plate data
% CPars.cuton (cuton frequency, in same units as s)
% CPars.phi  (phase angle of phase plate, in radians.  We then use max(alpha,phi)
% as the (negative) constant term in chi.
% 
% In the alternate form if the doComplex flag ==1, the complex result will
% be calculated, c=i*e^{i*pi*chi}
% 
% The frequency s can be of any dimension.
% Note that astigmatism makes sense only for 2D frequency values.
% For 1D or 3D frequency arrays, no astigmatism is applied, and the defocus
% is the value along the x-axis (ang=0) are taken.

deltadef=0;
cuton=0;
doComplex=0;
if isstruct(lambda)
    if nargin>2
        doComplex=defocus; % pick up third argument.
    end;
    P=lambda;
    lambda=P.lambda;
    defocus=P.defocus;
    Cs=P.Cs;
    B=P.B;
    alpha=P.alpha;
    if isfield(P,'cuton')
        cuton=P.cuton;
    end;
    if isfield(P,'phi')
        alpha=max(alpha,P.phi);
    end;
    if isfield(P,'deltadef')  % we are handling astigmatism
        deltadef=P.deltadef;
        theta=P.theta;
    end;
end;

if deltadef~=0 && ndims(s)<3  % handle astigmatism
    if all(size(s)>1)         % Truly 2-dimensional
        ang=atan2(s(2),s(1));
    else
        ang=0;
    end;
    defocus=defocus+deltadef*cos(2*(ang-theta));
end;
s2=s.^2;
chi=-1e4*lambda.*defocus.*s2+Cs.*lambda.^3.*5e6.*s2.^2-alpha/pi;

if doComplex
    c=1i*exp(1i*pi*chi-B*s2);
else
    c=sin(pi*chi).*exp(-B*s2);
end;
if cuton  % handle sharp cut-on of a phase plate.
    c=c.*(0.5+0.5*erf((abs(s)-cuton)*10/cuton));
end;
