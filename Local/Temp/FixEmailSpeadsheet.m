% FixEmailSpreadsheet.m
% generate a semicolon-delimited text file, by splitting off a trailing 
% email address in each line. This regenerated a NAMI support group address
% list.
fi=fopen('Emails.txt');
fo=fopen('EmailsConv.txt','w');
while ~feof(fi)
    li=fgetl(fi);
    s=textscan(li,'%s ');
    eAddr=s{1}{end};
    ptr=findstr(eAddr,li);
    eName=li(1:ptr-1);
    outLine=[eName ';' eAddr];
    disp(outLine);
    fprintf(fo,'%s\n',outLine);
end;
fclose(fi);
fclose(fo);    