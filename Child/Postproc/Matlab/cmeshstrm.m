function m=cmeshstrm(basenm,ts)
% CPLOTMESH: Plots Child triangulation.
%  filenm -- name of edge file
%  ts -- time slice to plot
%filesys=['/rhine/data/'];
%filesys=['/niagara/data/'];
filesys=['/arno/usr/users/'];
filenm= [filesys basenm '.nodes' ]
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.edges' ]
efid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ]
zfid=fopen(filenm,'r');
for i=1:ts
    tm = fscanf(nfid,'%f',1)
    tm = fscanf(efid,'%f',1);
    nn = fscanf(nfid,'%d',1)
    ne = fscanf(efid,'%d',1)
    e=fscanf(efid,'%f',[3,ne]); 
    n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
end
%hold off
e=e+1;
%plot3([0 0],[0 0],[0 0])
hold on
nx=rot90(n(1,:));
nx=rot90(nx);
nx=rot90(nx);
ny=rot90(n(2,:));
ny=rot90(ny);
ny=rot90(ny);
z=rot90(z);
z=rot90(z);
z=rot90(z);
m=[nx ny z];
filenm= [filesys basenm '.net' ]
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.area' ]
afid=fopen(filenm,'r');
for i=1:ts
  tm=fscanf(nfid,'%f',1) 
  tm=fscanf(afid,'%f',1);
  nn= fscanf(nfid,'%d',1)
  nn= fscanf(afid,'%d',1); 
  nbr=fscanf(nfid,'%d',[1,nn]);
  a=fscanf(afid,'%f',[1,nn]);
end
nd=m;
if nd==0, nd=creadxyz(basenm,ts); end
sz=size(nd);
if( sz(1)<nn ), error('Parameter ND must have at least NACTNODES columns');end
if( sz(2)~=3 ), error('Parameter ND must have 3 rows' );end
maxwid=4;
nbr=nbr+1;
lw=sqrt(sqrt(a));   % Line widths
%lw=sqrt(a);   % Line widths
lw=maxwid*lw/max(lw)+0.1;
hold on
%plot3([0 0],[0 0],[0 0])
amin = 50000;
%temporary: raise elevation so streams plot on top of mesh
%nd(:,3)=nd(:,3)+500;
for i=1:nn
  if a(i)>amin
    line( [nd(i,1) nd(nbr(i),1)], [nd(i,2) nd(nbr(i),2)], [nd(i,3) nd(nbr(i),3)],'color','k','linewidth',[lw(i)] ); 
    %set(h,'linewidth',[lw(i)]);
  end
end

% Now plot any sinks/lakes
%filenm= ['/rhine/data/gtucker/' basenm '.lakes' ]
%nfid=fopen(filenm,'r');
%if nfid>0
%  for i=1:ts
%    tm=fscanf(nfid,'%f',1)
%    nl= fscanf(nfid,'%d',1)
%    lk=fscanf(nfid,'%d',[1,nl]);
%  end
%  lk=lk+1;
%  for i=1:nl,plot3( [nd(lk(i),1)],  [nd(lk(i),2)],  [nd(lk(i),3)], 'g*' ),end
%end
%max(nd(:,3))
for i=1:2:ne
  plot3( [n(1,e(1,i)) n(1,e(2,i))], [n(2,e(1,i)) n(2,e(2,i))], [z(e(1,i)) z(e(2,i))],'k','linewidth',[0.1] ) 
  %line( [n(1,e(1,i)) n(1,e(2,i))], [n(2,e(1,i)) n(2,e(2,i))], 'color', 'y' ) 
end
hold off



