function sp=amSpectrum3(n,fs,amps)
ct=ceil((n+1)/2);
h=zeros(n,n);
f2=Radius(n).^2;
        for k=1:numel(fs)
            s=1/(2*pi*fs);
            a=2*pi*fs^2;
            h=sp+amps(k)*exp(-f2/fs(k)^2);
        end;
    end;
end;

