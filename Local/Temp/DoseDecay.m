% DoseDecay.m
% Compare Grigorieff and Rubinstein dose dependence
f=(.001:.001:.3)';
g=.75*(.245*f.^(-1.665)+2.81);
r=5+80./((f/.014).^2+1);

semilogy(f,[r g]);



