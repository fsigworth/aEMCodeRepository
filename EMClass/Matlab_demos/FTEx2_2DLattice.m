% FTEx2  FT of a lattice
% Click on the left panel to set the unit cell dimension.  The reciprocal
% lattice is shown on the right.  The keystroke q exits.
figure(1);
SetGrayscale;
n=256;
c=[n n]/2+1;  % center point
r=radius(n,c);
w=exp(-(2*r/n).^2);

g=zeros(16); g(9,9)=1;
g=GaussFilt(g,0.2);  % Gaussian function
p=ones(16);
m0=Mask(zeros(n),c, p,g);
subplot(1,2,1);
imacs(m0);
% mx=max(max(m0));
b=0;
while ~(b==32|b=='q') % space or 'y' exits the loop
    [dx, dy, b]=ginput(1);
    if b==1
        dx=max(round(abs(dx-c(1))),2);
        dy=max(round(abs(dy-c(2))),2);
        subplot(1,2,1);
        title([dx dy]); 
        subplot(1,2,2);
        title('...');
        drawnow;
        
        % make a lattice
        m=m0*0.5;  % Make the central point brighter than others.
        xr=floor(n/(2*dx))*dx;
        yr=floor(n/(2*dy))*dy;
        np=0;
        points=zeros(2,1);
        for x=c(1)-xr:dx:c(1)+xr
            for y=c(2)-yr:dy:c(2)+yr
                points(:,np+1)=[x y]';
                np=np+1;
            end;
        end;
        m=Mask(m,points,p,g);
        subplot(1,2,1);
        ms=w.*m;
        imacs(ms);
        title([dx dy]); drawnow;
        
        mf=fftshift(fftn(fftshift(ms))); 
        subplot(1,2,2);
        imacs(real(mf)); 
        title('FT'); drawnow;
    end;
end;

