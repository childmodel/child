function tri=ctrimesh(basenm,ts)
% CTRIMESH: Trimesh from Child triangulation
%  filenm -- name of edge file
%  ts -- time slice to plot
%         G. Tucker, 1998
filesys='';
filenm= [filesys basenm '.nodes' ]
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.tri' ]
tfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ]
zfid=fopen(filenm,'r');
for i=1:ts
  tm = fscanf(nfid,'%f',1)
  tm = fscanf(tfid,'%f',1);
  nn = fscanf(nfid,'%d',1)
  nt = fscanf(tfid,'%d',1)
  t=fscanf(tfid,'%f',[9,nt]); 
  n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
end
tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
x = n(1,:);
y = n(2,:);
trimesh(tri,x,y,z,zeros(nn,1))


