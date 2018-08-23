function str=EpochToDatestr(tUnix)
tMatlab = datenum (1970,1,1) + tUnix / 86400;
str=datestr(tMatlab);

