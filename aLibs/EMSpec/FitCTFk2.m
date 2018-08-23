function [P, cc, spect]=FitCTFk2(mc,Pa,pixA0,ds,fRange,options)
% function [P, cc, spect]=FitCTFk2(mc,Pa,pixA0,ds,fRange,options)
%  Derived from FitCTF, uses the algorithm of N. Grigorieff's CTFFIND3.
% With a proper options.carbonModel, it can compute an absolute B factor
% for a micrograph.  Fitted are all fields of Pa which have multiple
% elements, plus parameters const, amp and k2Dip which are returned as new
% fields in P.  k2Dip is the magnitude of the lf depression of the power
% spectrum in the K2 camera.
% Positive defocus values correspond to *underfocus*.
% Parameter ranges are defined by the structure Pa.  If vectors of values
% are given, all combinations of these values are tested in the brute-force
% search.  Afterwards, parameters are optimized using Simplex.
% The arguments pixA, maxres, minres are all in angstroms.
% The test argument controls display of results.  test=1 means graphical
% output (default).  test=2 gives text display of execution times.
% options is a struct of parameters for special options with fields
% spFromWholeMicrograph: use entire micrograph for spectrum.
%
% Here is an example:
%   Pa.lambda=EWavelength(200); % lambda in angstroms
%   Pa.defocus=0:0.5:20;  % in microns
%   Pa.deltadef=-2:.4:2;  % astigmatism
%   Pa.theta=[-90:30:90]*pi/180;
%   Pa.alpha=.02;  % This can also be fitted, for use with phase plate.
%       It is actually an angle in radians; to fit it, give a vector of
%       values e.g. 0.5:0.5:2.
%   Pa.B=20:20:100;  % This is fitted if a vector of values is given.
%
%    % The following parameters are not fitted
%   Pa.Cs=2;
%   Pa.lambda=.025
%    % Actual function call
%   P=FitCTFk2(m,Pa,1.4,2,[.05 .16],opts);
%   P      % print out the values

if nargin<6
    options=struct;
end;

% Set up default options

if ~isfield(options,'test')
    options.test=1;  % by default, show graphics.
end;
if ~isfield(options,'spFromWholeMicrograph')
    options.spFromWholeMicrograph=0;
end;
if ~isfield(options,'blockSize')
    options.blockSize=256;  % FFT block size, must be a multiple of 4.
end;
% if ~isfield(options,'fExponent')
%     options.fExponent=1;
% end;
if ~isfield(options,'k2Mode')
    options.k2Mode=1;
end;
if ~isfield(options,'carbonModel')
    options.carbonModel=ones(options.blockSize);
end;
if ~isfield(options,'title')
    options.title='';
end;

numRefs=2;
persistent R;  % We cache the references to allow faster execution of repeats.

mc=single(mc);
mc=mc-mean(mc(:));


pixA=pixA0*ds;  % pixel size in A.
Pa.pixA=pixA;

nu=options.blockSize;
f2d=RadiusNorm(nu)/ds;
df=1/(pixA*nu);  % spatial frequency unit of spectrum

% fExponent=options.fExponent;
k2Factor=sinc(f2d*2.58).^2;
k2Factor=exp(-(f2d/.231).^2.7);

% code for hp filtering of spectrum, presently unused.
% defocusMin=min(Pa.defocus);
% firstZero=sqrt(defocusMin*Pa.lambda*1e4); % first zero, in Å
% w0=firstZero/2;  % High-pass filter for spectrum, in A.  Should be smaller than
% w0=100;
% w=w0/pixA;     % width of real-space kernel.

freqs=(0:df:df*(nu/2-1))';

% ---------- Compute the power spectrum -------------
[nx, ny]=size(mc);
nv=round(nu/4);   % overlap
tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction

if options.test>0  % Show the micrograph
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
    imacs((1:nx)*pixA/10,(1:ny)*pixA/10,BinImage(mc,4));
    xlabel('nm');
    title(options.title,'interpreter','none');
    drawnow;
end;

sd=zeros(tx,ty);
% Make a window for spectral estiation
window=fuzzymask(nu,2,0.45*nu,0.1*nu);
winsum=sum(window(:));

sps=zeros(nu,nu,tx,ty);

if options.test>1
    disp('Computing spectra...');
    tic
end;
for ix=1:tx
    [x0 x1 u1]=TileCoords(nx,nu,nv,ix);
    
    for iy=1:ty
        [y0 y1 v1]=TileCoords(ny,nu,nv,iy);
        % I don't think it's necessary to remove gradients, but Niko does
        % this:
        tm=RemoveGradients(double(mc(x0+1:x0+nu,y0+1:y0+nu)));
        tm=tm.*window;
        tm=tm-sum(tm(:))*window/winsum;
        
        sp2=abs(fftn(tm)).^2;
        sps(:,:,ix,iy)=sp2;
        sd(ix,iy)=sqrt(sum(sp2(:)))/nu;
    end;
end;
if options.test>1
    toc
end;

sdmn=min(sd(:));
sdmx=max(sd(:));
if options.spFromWholeMicrograph
    sdmin=sdmn;
    sdmax=sdmx;
else
    % Make a histogram of the s.d. values and find tiles having s.d. near the
    % mode.
    nbins=tx+2;
    hthresh=0.2;
    % hthresh=-1;  % no selection
    
    dbin=(sdmx-sdmn)/(nbins-1);
    bins=sdmn:dbin:sdmx;
    
    [h, x]=histc(sd(:),bins);
    
    [mxh, j]=max(h);
    % j is the mode of the distribution
    % Find the first bin below hthresh on the high side
    jmax=j-1+find(h(j:nbins)<hthresh*mxh,1,'first');
    if numel(jmax)<1
        jmax=nbins;
    end;
    
    % find the bin below hthresh on the low side
    jmin=find(h(1:j)<hthresh*mxh,1,'last');
    if numel(jmin)<1
        jmin=1;
    end;
    sdmin=bins(jmin);
    sdmax=bins(jmax)+dbin;
end;
% Show histogram
% subplot(2,3,2);
% bar(x,h);
% hold on
% plot(sdmin,0,'w.');
% plot(sdmax,0,'w.');
% hold off

% Show regions that we used.
if options.test
    subplot(2,3,2);
    imacs((sd<=sdmax).*(sd>=sdmin));
    drawnow;
end;

cumsp=zeros(nu,nu);
count=0;
for ix=1:tx
    for iy=1:ty
        if (sd(ix,iy)<=sdmax) && (sd(ix,iy)>=sdmin)
            cumsp=cumsp+sps(:,:,ix,iy);
            count=count+1;
        end;
    end;
end;
if options.test>1
    disp([num2str(count) ' tiles of ' num2str(tx*ty) ' used']);
end;
spect=fftshift(single(cumsp/count));

%% ---------- mask and normalize the spectrum -----------

sp1=fftshift(cumsp);
% determine resolution limits
rmin=fRange(1)/df;
dmin=round(rmin);
rmax=fRange(2)/df;
dmax=min(nu/2,round(rmax));

outerlimits=fuzzymask(nu,2,rmax,rmax/20);
limits=outerlimits-fuzzymask(nu,2,rmin,1);  % annulus from rmin to rmax

% get the mean of the high-frequency end of the 2d spectrum
annulus=fuzzymask(nu,2,0.4*nu,0.1*nu)-fuzzymask(nu,2,0.3*nu);  % from 60 to 80% of Nyquist

spMean=(annulus(:)'*sp1(:))/sum(annulus(:));
spectN=sp1/spMean;  % normalized to have HP part approx =1.
spectMsk=limits.*sp1/spMean;  % same, but with limits imposed.

halfLimits=limits(nu/2+1:nu,:);  % right half-plane
halfLimits(1,nu/2:nu)=0;  % zero out the upper meridian
halfSp=spectMsk(nu/2+1:nu,:);
halfSpv=halfSp(:);



%% -------------- Get references ----------------
% % Here we test for an existing set of references, and compute new ones
% % only if necessary.
% First, check whether the static variable R exists at all.
b=0;
if exist('R','var')
    b=numel(R);
end;

% Now see if the problem has changed from previous calls.
ind=0;
for i=1:b
    if StructsEqual(R(i).Pa0,Pa) && size(R(i).refs,1)==numel(halfSpv)
        ind=i;
    end;
end;

if ind==0
    disp('Making new references');
    ind=min(b+1,numRefs);
    if options.test>1
        disp('Making new CTF references');
    end;
    
    [R(ind).refs, R(ind).refsz]=MakeCTFRefs(nu,pixA,Pa,halfLimits);
    %     [R(ind).refso R(ind).refsz]=MakeCTFRefs(nu,res,Pao,halflimits);
    
    R(ind).refs=reshape(R(ind).refs,nu^2/2,prod(R(ind).refsz));
    %     R(ind).refso=reshape(R(ind).refso,nu^2/2,prod(R(ind).refsz));
    R(ind).Pa0=Pa;
    if options.test>1
        whos R
    end;
end;

%% ----------- Initial fit --------------

if options.test>1
    disp('Cross-correlation')
end;

nRefs=prod(R(ind).refsz);
ccx=zeros(nRefs,1);
lsPars=zeros(3,nRefs);
F=zeros(nu^2/2,3);
F(:,1)=halfLimits(:);  % constant
q=k2Factor(nu/2+1:nu,:);  % k2Factor term
F(:,2)=reshape(q.*halfLimits,nu^2/2,1);
halfCModel=options.carbonModel(nu/2+1:nu,:);
halfCModel=halfCModel(:);
for ir=1:nRefs
    F(:,3)=R(ind).refs(:,ir).*halfCModel;
    a=LinLeastSquares(F,halfSpv);
    model=F*a;
    diff=halfSpv-model;
    ccx(ir)=diff(:)'*diff(:);
    lsPars(:,ir)=a;
end;
% % nrefs=size(R(ind).refs,2);
% cc=R(ind).refs'*halfspv;  % Cross correlation done here!
%
% pick up the indices of the best match
[mxc, mxi]=max(-ccx);
a=lsPars(:,mxi);

ccx=reshape(-ccx,R(ind).refsz);
[m, l, i, j, k]=ind2sub(R(ind).refsz,mxi);

if options.test>1
    Initial_c=mxc
end;

P=Pa;
% Set the starting search values
P.defocus=Pa.defocus(i);
P.deltadef=Pa.deltadef(j);
P.theta=Pa.theta(k);
P.alpha=Pa.alpha(l);
P.B=P.B(m);
P.const=a(1);
P.k2Dip=a(2);
P.amp=a(3);

if options.test
    subplot(2,3,3);
% % % % % %     ContourPlot(Pa,ccx,i,j,k,l,m);
    % Make the split-power spectrum display
    SplitDisplays(spectN,P,pixA,1);
    drawnow;
end;

%% Do the optimization
% disp('Final optimization');
if options.test>1
    disp('Final optimization.');
end;

% Set up the fields for parameters to vary
fields=fieldnames(Pa);
nf=numel(fields);
aFields={};
aPars=[];
msk=P;
naf=0;
extraFields={'const','k2Dip','amp'};
for i0=1:nf
    fname=fields{i0};
    var=numel(Pa.(fname))>1;
    msk.(fname)=var;
    if var
        naf=naf+1;  % number of active fields
        aFields{naf}=fname;  % active fields
        aval=P.(fname);
        aPars(naf)=aval; % active parameters
    end;
end;
for i0=1:numel(extraFields)
    
    naf=naf+1;
    fname=extraFields{i0};
    aFields{naf}=fname;
    aPars(naf)=P.(fname);
end;

% Pinit=P

opts=optimset('TolX',1e-6,'FinDiffRelStep',1e-3,...
    'display','off');

[pars,cc,residual,exitFlag,output]=lsqnonlin(@EvalNCC,aPars,[],[],opts);

% Put the parameters back into a structure
for i0=1:naf
    P.(aFields{i0})=pars(i0);
end;

if options.test
    SplitDisplays(spectN,P,pixA,1);
    title(['CC =' num2str(cc)]);
end;


if options.test>1
    cc
    Final_values=[P.defocus P.deltadef 180/pi*P.theta P.alpha]
end;  % end of main function





    function y=EvalNCC(x)
        %         disp(x(:)')
        Pt=P;
        for ifld=1:naf
            Pt.(aFields{ifld})=x(ifld);
        end;
        fcn=limits.*EvalModel(Pt);
        y=double(spectMsk(:)-fcn(:));
    end


    function fcn=EvalModel(Pt)
        ct2=CTF(nu,pixA,Pt).^2.*limits.*options.carbonModel;
        %         fcn=Pt.amp*ct2+Pt.const*(1-Pt.k2Dip*exp(-(f2d/.231).^2.7));
        fcn=Pt.amp*ct2+Pt.const+Pt.k2Dip*k2Factor;
    end


    function SplitDisplays(sp,P,pixA,showPars)
        % Make the normalized spectrum and fit
        subplot(2,3,6);
        Ptm=P;
        %     Ptm.k2Dip=0;
        ref2=EvalModel(Ptm);
        ref1=RadialCTF(ref2,P,pixA);
        
        %     Ptm=P;
        %     Ptm.amp=0;
        %     dipModel=EvalModel(Ptm);
        %     dipModel=1;
        
        %     spAnn=annulus.*sp;
        %     spHFMean=sum(spAnn(:))/sum(annulus(:));
        %     nsp2=sp./(dipModel*spHFMean);
        nsp2=sp;
        nsp1=RadialCTF(nsp2,P,pixA);
        
        %     Make the 1d plot
        specs1D=[nsp1 ref1];
        plot(freqs(dmin:dmax),specs1D(dmin:dmax,:));
        xlabel('A^{-1}');
        ylabel('Spectral density');
        if showPars
            title(['d = ' num2str(P.defocus,4) '  \Deltad = ' num2str(P.deltadef,3)...
                '  \Theta = ' num2str(180/pi*P.theta,3) '  B = ' num2str(P.B,3)]);
        else
            title('...');
        end;
        %     Make 2d displays
        limV=limits(:);
        meRef=(limV'*ref2(:))/(limV'*limV);
        meNsp=(limV'*nsp2(:))/(limV'*limV);
        me=(meRef+meNsp)/2;
        mePad=(1-limits)*me;
        ref2fl=limits.*ref2+mePad;
        nsp2fl=limits.*nsp2+mePad;
        img2=nsp2fl;
        img2(1:nu/2+1,:)=ref2fl(1:nu/2+1,:);  % create the split image
        
        subplot(2,3,5);
        xs=(-nu/4:nu/4-1)*df;
        imacs(xs,xs,Crop(img2,nu/2));
        xlabel('A^{-1}');
        title(['\alpha = ' num2str(P.alpha*180/pi,4) ' deg']);
        
        subplot(2,3,4);
        xs=(-nu/2:nu/2-1)*df;
        imacs(xs,xs,img2);
        xlabel('A^{-1}');
        
    end

    function ContourPlot(Pa,ccx,i,j,k,l,m)
        % make a contour plot of two of the parameters.
        % inidices of refsz are(nB nalpha ndef ndeldef ntheta)
        %     We don't plot theta, but use all the others
        w=shiftdim(ccx(:,:,:,:,k),1);  % w(alpha, def, deldef, theta, B)
        inds=[l i j m];  % peak indices of w
        fields={'alpha' 'defocus' 'deltadef' 'B'};
        % pick the first two to plot
        sz=zeros(numel(fields),1);
        for id=1:numel(fields)
            sz(id)=numel(Pa.(fields{id}));
        end;
        ip=find(sz>1);  % find the first two fields containing vectors
        while numel(ip)>2  % freeze all later indices
            p1=ip(end)-1;
            ws=(shiftdim(w,p1));
            w=shiftdim(ws(inds(p1+1),:,:,:),4-p1);
            ip(end)=[];
        end;
        w=squeeze(w);
        
        switch numel(ip)  % draw nothing is nothing is variable.
            case 1  % only one variable parameter, make a plot
                plot(Pa.(fields{ip}),w);
                xlabel(fields{ip});
                ylabel('Correlation coefficient');
                drawnow;
            case 2
                contourf(Pa.(fields{ip(2)}),Pa.(fields{ip(1)}),w,10);
                xlabel(fields{ip(2)});
                ylabel(fields{ip(1)});
                %                 colorbar;
                drawnow;
        end;
        
    end
end