function [x,h,xp,yp,sp,ap]=cstrmproseries( basenm, ts, x, y )
% CSTRMPROSERIES: Retrieves and plots a series of longitudinal stream 
%               profiles from a CHILD run, starting from a given node at
%               location (x,y), and moving downstream from there until an
%               open boundary is reached. Returns vectors for the height and
%               distance along the profile. Same as cstreamprofilexy but
%               plots multiple profiles through time.
%
%  Usage: [dist,height,xp,yp,sp,ap] = cstrmproseries( basenm, ts, x, y )
%
% Parameters:
%  basenm -- name of files with run
%  ts -- time slice to plot
%  x,y -- coordinates of starting point
%
% Returns: downstream distance, height, coordinates of each point, and
%          slope and drainage area at each point
%
% Functions called: creadxyzb.m, cplotmesh.m, cplotarea.m, cread.m,
%                   cfindstrmprofile.m
%
%   G. Tucker, 2002; 2008
%

% INITIALIZE
% Open files and handle file-access errors
dirfid = fopen([basenm '.net'],'r');
if dirfid<=0,error(['Unable to open flow-direction file ' basenm '.net']);end
nfid=fopen([basenm '.nodes'],'r');
if dirfid<=0
    fclose(dirfid);
    error(['Unable to open node file ' basenm '.nodes']);
end
zfid=fopen([basenm '.z'],'r');
if zfid<=0
    fclose(dirfid);
    fclose(nfid);
    error(['Unable to open elevation file ' basenm '.z']);
end
sfid=fopen([basenm '.slp'],'r');
if sfid<=0
    fclose(zfid);
    fclose(dirfid);
    fclose(nfid);
    error(['Unable to open slope file ' basenm '.slp']);
end
afid=fopen([basenm '.area'],'r');
if afid<=0
    fclose(sfid);
    fclose(zfid);
    fclose(dirfid);
    fclose(nfid);
    error(['Unable to open drainage-area file ' basenm '.area']);
end



% PROCESS
% Cycle through time slices, reading time-slice data from each file and
% sending the data to cstrmprofilexy to generate a profile
for i=1:ts
  
  % Read x,y,e,b data from node file
  tm = fscanf(nfid,'%f',1);
  nn = fscanf(nfid,'%d',1);
  n=fscanf(nfid,'%f',[4,nn]);
  
  % Read z data from elevation file
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
  
  % Read flow-direction data from direction file
  tm = fscanf(dirfid,'%f',1);
  nn = fscanf(dirfid,'%d',1); 
  dir=fscanf(dirfid,'%f',[1,nn]);
  dir=dir+1;  % Change of indices from 0..N-1 to 1..N
  
  % Read slope data from slope file
  tm = fscanf(sfid,'%f',1);
  nn = fscanf(sfid,'%d',1); 
  s=fscanf(sfid,'%f',[1,nn]);
  
  % Read drainage area data from area file
  tm = fscanf(afid,'%f',1);
  nna = fscanf(afid,'%d',1); 
  atmp=fscanf(afid,'%f',[1,nna]);
  a=0*s;
  a(1:nna)=atmp;
  atmp=[];

  % Configure xyzb data
  xyzb = [ rot90(n(1,:),3) rot90(n(2,:),3) rot90(z,3) rot90(n(4,:),3) ];

  distmax = 0.0;  % Maximum streamwise distance found so far

  dx=xyzb(:,1)-x;
  dx=dx.*dx;
  dy=xyzb(:,2)-y;
  dy=dy.*dy;
  dist=sqrt(dx+dy);
  [mindist index]=min(dist);
  %figure(2)
  [dist,h,xp,yp,sp,ap]=cfindstrmprofile( xyzb, dir, index, s, a );
  if max(dist) > distmax, distmax = max(dist); end
  plot( dist, h )
  hold on

fprintf('Maximum profile distance is %f\n',distmax);
  
end






