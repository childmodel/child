function [tri,meanElev]=ctrisurf(basenm,ts)
% CTRISURFTX: Trisurf from Child triangulation. Colored by surface sediment
%             texture (%sand)
%  filenm -- name of output file
%  ts -- time slice to plot
%    Usage: [t,z]=ctrisurf(filenm,ts), where t=triangles, z=mean elevations
%    G. Tucker, 2000 
filesys=['/arno/usr/users/gtucker/'];
%filesys=['/rhine/data/'];
%filesys=['/niagara/usr/users/gtucker/'];
%filesys=['/niagara/data/'];
filenm= [filesys basenm '.nodes' ]
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.tri' ]
tfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ]
zfid=fopen(filenm,'r');
filenm= [filesys basenm '.tx' ];
txfid = fopen(filenm,'r');
if txfid<=0, error(['Cannot open ' filenm]);end
meanElev=[];
for i=1:ts
  tm = fscanf(nfid,'%f',1)
  tm = fscanf(tfid,'%f',1);
  tm = fscanf(txfid,'%f',1);
  nn = fscanf(txfid,'%d',1);
  nn = fscanf(nfid,'%d',1);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 
  n=fscanf(nfid,'%f',[4,nn]);
  tx=fscanf(txfid,'%f',[1,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
  meanElev = [meanElev mean(z)]; 
end
tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
x = n(1,:);
y = n(2,:);
trisurf(tri,x,y,z,tx)
rot90(meanElev,3);

