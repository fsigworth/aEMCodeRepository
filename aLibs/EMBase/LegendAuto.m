function LegendAuto(n)
%function LegendAuto(n)
% Make a legend with n integer entries.  For example
% plot(m);
% LegendAuto(size(m,2));
% puts up a legend entry for each column of m.
ns=1:n;
legend(num2str(ns'));

