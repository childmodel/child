function a = chsa(filenm,time)
%CHSA: CHILD slope-area plot
%    Usage: chsa( file name, time slice # )
%       G. Tucker, 1997
filenm= ['/nimbus/staff/academic/gtucker/Dev/Child/Examples/' filenm]
sfid = fopen([filenm '.slp'],'r');
afid = fopen([filenm '.area'],'r');
for i=1:time
  tm=fscanf(sfid,'%f',1)
  np=fscanf(sfid,'%d',1);
  tm=fscanf(afid,'%f',1);
  np=fscanf(afid,'%d',1);
  s=fscanf(sfid,'%f',[np]);
  a=fscanf(afid,'%f',[np]);
end
loglog(a,s,'.')
xlabel('Drainage area (m2)')
ylabel('Slope')
title(filenm)


