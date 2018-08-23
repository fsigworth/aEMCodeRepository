% TextFileToCellArray.m

[name,pa]=uigetfile('*.txt');
cd(pa);
f=fopen(name);
q=textscan(f,'%s');
allNames=q{1};
