function m=cplotarea(basenm,nd,ts,amin)
% CPLOTAREA: Plots network from Child simulation, with line thickness in 
% proportion to drainage area. Call this after running CPLOTMESH.
%  Usage: cplotareard( filename, nd, ts, {amin} )
% Parameter ND is an NNODES by 3 array containing the x,y, and z coords of
% the nodes. It's obtained from CPLOTMESH.
%      G. Tucker, 1997
filesys=[''];
filenm= [filesys basenm '.net' ];
nfid=fopen(filenm,'r');
if nfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.area' ];
afid=fopen(filenm,'r');
%nd(:,3)=nd(:,3)+5;
fprintf('CPLOTAREA: Reading data ...\n');
if nd==0
    nd=creadxyz(basenm,ts);
end
if nargin<4
    amin=0;
end
for i=1:ts
  tm=fscanf(nfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(afid,'%f',1);
  nn= fscanf(nfid,'%d',1);
  nn= fscanf(afid,'%d',1); 
  nbr=fscanf(nfid,'%d',[1,nn]);
  a=fscanf(afid,'%f',[1,nn]);
end
fprintf('\n');
fprintf('CPLOTAREA: Plotting...');
sz=size(nd);
if( sz(1)<nn ), error('Parameter ND must have at least NACTNODES columns');end
if( sz(2)~=3 ), error('Parameter ND must have 3 rows' );end
maxwid=5;
nbr=nbr+1;
%lw=sqrt(sqrt(a));   % Line widths
lw=sqrt(a);   % Line widths
%lw=a;
lw=maxwid*lw/max(lw)+1;
plot3([0 0],[0 0],[0 0])
hold on
for i=1:nn 
  if i==nbr(i),i,nbr(i),end
  if a(i)>amin && nd(i,3)>0
    h=plot3( [nd(i,1) nd(nbr(i),1)], [nd(i,2) nd(nbr(i),2)], [nd(i,3) nd(nbr(i),3)]+2,'k' ); 
    set(h,'linewidth',[lw(i)]);
  end
end
hold off

% Now plot any sinks/lakes
% filenm= ['/rhine/data/gtucker/' basenm '.lakes' ];
% nfid=fopen(filenm,'r');
% if nfid>0
%   for i=1:ts
%     tm=fscanf(nfid,'%f',1);
%     nl= fscanf(nfid,'%d',1);
%     lk=fscanf(nfid,'%d',[1,nl]);
%   end
%   lk=lk+1;
%   for i=1:nl,plot3( [nd(lk(i),1)],  [nd(lk(i),2)],  [nd(lk(i),3)], 'g*' ),end
% end
% 
m=nbr; %debug
fprintf('\n');

