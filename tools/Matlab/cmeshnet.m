function m = cmeshnet( filename, ts )
%CMESHNET
m = cplotmeshm2( filename, ts )
cplotaream2( filename, m, ts )
view(0,90)
set(gca,'dataaspectratio',[1 1 1])
