function fm=sft(m);
fm=fftshift(fftn(fftshift(m)));
