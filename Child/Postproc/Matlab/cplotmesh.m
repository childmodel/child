function m=cplotmesh(basenm,ts)
% CPLOTMESH: Plots Child triangulation.
%  filenm -- name of edge file
%  ts -- time slice to plot
filesys='';
filenm= [filesys basenm '.nodes' ]
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.edges' ]
efid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ]
zfid=fopen(filenm,'r');
fprintf('CPLOTMESH: Reading data...\n');
for i=1:ts
    tm = fscanf(nfid,'%f',1);
    fprintf('Time slice %d (T=%f)\n',i,tm);
    tm = fscanf(efid,'%f',1);
    nn = fscanf(nfid,'%d',1);
    ne = fscanf(efid,'%d',1);
    e=fscanf(efid,'%f',[3,ne]); 
    n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
end
fprintf('\n');
hold off
e=e+1;
plot3([0 0],[0 0],[0 0])
hold on
for i=1:2:ne
   h=plot3( [n(1,e(1,i)) n(1,e(2,i))], [n(2,e(1,i)) n(2,e(2,i))], [z(e(1,i)) z(e(2,i))], 'k' );
  %line( [n(1,e(1,i)) n(1,e(2,i))], [n(2,e(1,i)) n(2,e(2,i))], 'color', 'k' ) 
%     set(h(1),'linewidth',[0.1])
end
drawnow
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


