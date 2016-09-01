function [tri,meanElev]=ctrisurf(basenm,ts,c,optfixedmesh)
% CTRISURF: Trisurf from Child triangulation
%  filenm -- name of output file
%  ts -- time slice to plot
%  c -- field to color (optional; default = color by elevation)
%    Usage: [t,z]=ctrisurf(filenm,ts,{c},{optfixedmesh}), where t=triangles, z=mean 
%                   elevations
%    Note: if you want to color by something other than elevation, read 
%         that data field in (for the same time step) using CREAD and give 
%         it in parameter 'c'. 
%    G. Tucker, 1998; 2008

%a=rand(3)
if nargin<4
    optfixedmesh = 0;
end
filesys='';
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');
if nfid<=0, error('Unable to open node file'),end
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0, fclose(nfid); error('Unable to open triangle file'),end
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
if zfid<=0, fclose(nfid); fclose(tfid); error('Unable to open elevation file'),end
meanElev=[];
for i=1:ts
  tm = fscanf(zfid,'%f',1);
  fprintf('CTRISURF: Reading time %f\n',tm);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
  if optfixedmesh==0 || i==1
      tm = fscanf(nfid,'%f',1);
      tm = fscanf(tfid,'%f',1);
      nn = fscanf(nfid,'%d',1);
      nt = fscanf(tfid,'%d',1);
      t=fscanf(tfid,'%f',[9,nt]);
      n=fscanf(nfid,'%f',[4,nn]);
  end
  meanElev = [meanElev mean(z)]; 
  if feof(nfid) || feof(tfid) || feof(zfid)
     fclose(nfid);
     fclose(tfid);
     fclose(zfid);
     error(['Reached end of file: there may be only ' num2str(i-1) ' output steps.']);
  end
end
fclose(nfid);
fclose(tfid);
fclose(zfid);

tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
x = n(1,:);
y = n(2,:);
b = n(4,:);
if nargin<3 %|| length(c)<length(z)
    c=z; 
elseif c==0
    c=z;
end
c = cinterpedges( x, y, c, b );  % this removes edge effects for plotting
if size(z,2)~=nn, error('Color index (c) must have same size as # nodes'),end
trisurf(tri,x,y,z,c)
rot90(meanElev,3);
