function m=creadxyzb(basenm,ts)
% CREADXYZ: Reads node x,y,z coordinates and boundary flag ("b")
%           for one timestep of a CHILD run
%           and returns them in a Nx4 matrix, where N=number of nodes.
%
%  Usage: m = creadxyzb( basenm, ts )
%
%  Parameters:
%    filenm -- name of edge file
%    ts -- time slice # to fetch data for
%
filesys=[''];
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');
if nfid<1
    error(['In creadxyzb.m, unable to open ' basenm '.nodes'])
end
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
for i=1:ts
  tm = fscanf(nfid,'%f',1);
  nn = fscanf(nfid,'%d',1);
  n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
end
m = [ rot90(n(1,:),3) rot90(n(2,:),3) rot90(z,3) rot90(n(4,:),3) ];
fclose(nfid);
fclose(zfid);


