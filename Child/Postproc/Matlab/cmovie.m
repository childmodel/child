function m=cmovie(basenm,nframes,ax,ay,az,ac)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%    Usage: m = cmovie( basenm, nframes, ax, ay, az, ac )
%
% Parameters:
%   basenm - base file name of CHILD run
%   nframes - number of frames
%   ax (optional) - extent of x-axis
%   ay (optional) - extent of y-axis
%   az (optional) - extent of z-axis
%   ac (optional) - upper end of color axis
%

% Open files
filesys='';
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');
if nfid<=0, error('Unable to open node file'),end
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0, error('Unable to open triangle file'),end
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
if zfid<=0, error('Unable to open elevation file'),end

for i=1:nframes
  
  % Read stuff
  tm = fscanf(nfid,'%f',1);
  fprintf('CTRISURF: Reading time %f\n',tm);
  tm = fscanf(tfid,'%f',1);
  nn = fscanf(nfid,'%d',1);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 
  n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);  
  
  % Extract coordinates
  tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
  x = n(1,:);
  y = n(2,:);
  b = n(4,:);

  % Remove edge effects
  %z = cinterpedges( x, y, z, b );  % this removes edge effects for plotting
  
  % Plot
  trisurf(tri,x,y,z,z);
  
  % Set plotting stuff
  if nargin>=5
    axis([0 ax 0 ay 0 az])
  end
  if nargin>=6
    caxis([0 ac])
  end
  %view(-30,60)
  %shading interp
  %colormap pink
  %colorit(2)
  
  % Capture the frame
  m(i)=getframe;
  
end

% Play the move
movie(m);

% close the files
fclose(nfid);
fclose(tfid);
fclose(zfid);
