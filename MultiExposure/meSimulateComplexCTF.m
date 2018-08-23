function filteredM = meSimulateComplexCTF(m,mi,ds,weights)
% % Given an object m (units are phase shift in radians) compute the output
% % of the image-merging process according to the mi.ctf, with weights being
% % the weights given to the respective images, by default [1 1 1...].
% % Note that this gives a signal corresponding to 2*CTF computed
% classically.
% 
n=size(m);
if nargin<3 || ds==0
    ds=mi.imageSize(1)/n(1);
end;
if nargin<4
    weights=mi.doses>0;
end;

effPixA=mi.pixA*ds;
freqs=RadiusNorm(n)/effPixA;  % frequencies for evaluating the CTF

% Compute the decay constant due to radiation, as a function of frequency
% --a fit to the data in the Rubenstein paper.
n0=2*0.16./(abs(freqs)+.004)+5;  % twice the critical dose at 200 kV
doses=single(mi.doses(:));
cumdose=cumsum([0; doses]);
q=find(doses,1);
dose1=doses(q);

coeffs=meComputeMergeCoeffs(freqs, mi.ctf, mi.doses,1,weights);

fTotal=complex(zeros(n));
indices=find(weights(:)'>0);
for i=indices
    P=mi.ctf(i);
    P.B=0;  % we handle the B factor after computing the intensity
    signal=(exp(-cumdose(i)./n0).*n0.*(1-exp(-doses(i)./n0)))/doses(i);
%     Assume that radiation damage blurs the specimen itself.  The
%     alternative model would be that radiation damage decoheres the
%     wavefront.
    decayedM=real(ifftn(fftn(m).*ifftshift(signal)));
    fm=fftn(exp(1i*exp(-1i*P.alpha)*decayedM));  % complex wavefront
    ctfx=ComplexCTF(n,effPixA,P);
    I=abs(ifftn(fm.*ifftshift(ctfx))).^2-1;  % Intensity of wavefront
    env=exp(-freqs.^2.*mi.ctf(i).B);  % B factor and CCD
    fTotal=fTotal+fftn(I).*ifftshift(env.*coeffs(:,:,i))*doses(i)/dose1;
end;
filteredM=real(ifftn(fTotal.*ifftshift(CCDSqrtDQE(n,ds))));
