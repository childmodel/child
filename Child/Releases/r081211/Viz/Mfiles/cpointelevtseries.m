function [t,zp]=cpointelevtseries(basenm,xp,yp,ntimes,amin)
% CPOINTELEVTSERIES: extracts a time series of elevation at/near a
% specified point. If desired, the closest point with a drainage area above
% a given threshold will be found
%    Usage: [t,z] = cpointelevtseries( basenm, xp, yp, {ntimes}, {amin} )
%
% GT, July 2008

if nargin<4, ntimes=101; end
if nargin<5, amin=0; end

% Open files
filenm= [basenm '.nodes' ];
nfid=fopen(filenm,'r');
if nfid<=0, error('Unable to open node file'),end
filenm= [ basenm '.z' ];
zfid=fopen(filenm,'r');
if zfid<=0, error('Unable to open elevation file'),end

if amin>0
    filenm= [ basenm '.area' ];
    afid=fopen(filenm,'r');
    if afid<=0, error('Unable to open area file'),end
end

cur_time = 1;

t=[];
zp=[];

while cur_time<=ntimes && ~feof(nfid)
  
  % Read stuff
  tm = fscanf(nfid,'%f',1);
  fprintf('CPOINTELEVTSERIES: Reading time %f\n',tm);
  nn = fscanf(nfid,'%d',1);
  n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);  
  
  t=[t tm];

  % Find the closest node to point
  x=n(1,:);
  y=n(2,:);
  if amin>0
      tm = fscanf(afid,'%f',1);
      na = fscanf(afid,'%d',1);
      a=fscanf(afid,'%f',[1,na]);
      x=x(a>amin);
      y=y(a>amin);
      z=z(a>amin);
  end
  delx = x-xp;
  dely = y-yp;
  dist_sq = delx.^2+dely.^2;
  [minval,closest_node]=min(dist_sq);
  fprintf('Time %f: Closest node is %f at (%f,%f), a distance %f away.\n',cur_time,closest_node,x(closest_node),y(closest_node),sqrt(minval));
  zp=[zp z(closest_node)];
  
  cur_time=cur_time+1;
  
end

% close the files
fclose(nfid);
fclose(zfid);
fclose(afid);
