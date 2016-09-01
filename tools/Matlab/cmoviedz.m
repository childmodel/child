function mdz = cmoviedz(filenm,time)
%CMOVIEDZ: Not really a movie per se, but a sequence of frames displayed
%          showing topography colored by depth of erosion/deposition since
%          previous time slice in CHILD run.
%  Usage:  mdz = cmoviedz( filename, number of frames )
%     G. Tucker, 1998
%filesys='/rhine/data/';
filesys='/niagara/data/';
fullnm= [filesys filenm '.nodes']
xyfid = fopen(fullnm,'r');
fullnm= [filesys filenm '.z']
zfid = fopen(fullnm,'r');
fullnm= [filesys filenm '.tri' ]
tfid=fopen(fullnm,'r');
colormap jet
flipcolormap
mindzv=[];
maxdzv=[];
meandzv=[];
% Read first time slice
tm0=fscanf(zfid,'%f',1);
tm=fscanf(xyfid,'%f',1);
tm = fscanf(tfid,'%f',1);
nt = fscanf(tfid,'%d',1);
nn=fscanf(xyfid,'%d',1);
nnz=fscanf(zfid,'%d',1);
z=fscanf(zfid,'%f',[1,nnz]);
xy=fscanf(xyfid,'%f',[4,nn]);
t=fscanf(tfid,'%f',[9,nt]); 
% Read remaining time slices
for i=2:time
  tm=fscanf(xyfid,'%f',1);
  tm=fscanf(zfid,'%f',1)
  tm = fscanf(tfid,'%f',1);
  dtm = tm-tm0
  tm0=tm;
  nt = fscanf(tfid,'%d',1);
  nn=fscanf(xyfid,'%d',1);
  nnz=fscanf(zfid,'%d',1);
  if nn~=nnz, error('Different no. of nodes in the two files.'),end
  xy0=xy;
  xy=fscanf(xyfid,'%f',[4,nn]);
  t=fscanf(tfid,'%f',[9,nt]); 
  z0=z;
  z=fscanf(zfid,'%f',[1,nn]);
  dz=getcdz(xy,xy0,z,z0);
  dz=dz/dtm;
  % correct for uplift:
  uplift=0;
  %uplift=0.00002*dtm;
  dz=dz-uplift;
  mindz=min(dz);
  maxdz=max(dz);
  meandz=mean(dz);
  %For display:
  %scalemax=0.005;   % Max value for scaling colors, for doing a series
  scalemax=0.1793;
  mxx=max([abs(mindz) maxdz scalemax] );
  dz(1)=mxx;
  dz(2)=-mxx;
  tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
  x = xy(1,:);
  y = xy(2,:);
  trisurf(tri,x,y,z,dz)
  view(60,80)
  %view(0,90)
  %axis([0 7000 0 7000 0 200])
  colorbar
  title(num2str(tm))
  mindz
  maxdz
  meandz
  mindzv=[mindzv mindz];
  maxdzv=[maxdzv maxdz];
  meandzv=[meandzv meandz];
  pause
end
min(mindzv)
mdz=[rot90(mindzv,3) rot90(maxdzv,3) rot90(meandzv,3)];

