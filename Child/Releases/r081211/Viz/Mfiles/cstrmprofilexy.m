function [x,h,xp,yp,sp,ap]=cstrmprofilexy( basenm, ts, x, y )
% CSTRMPROFILEXY: Retrieves and plots a longitudinal stream profile from a
%               CHILD run, starting from a given node which the user selects
%               with the mouse, and moving downstream from there until an
%               open boundary is reached. Returns vectors for the height and
%               distance along the profile. Same as cstreamprofile but
%               takes starting coordinates as inputs
%
%  Usage: [dist,height,xp,yp,sp,ap] = cstrmprofile( basenm, ts, x, y )
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

dir=cread([basenm '.net'],ts);
dir=dir+1;  % Change of indices from 0..N-1 to 1..N
xyzb = creadxyzb(basenm,ts);

% slope and area
slp = cread([basenm '.slp'],ts);
artmp = cread([basenm '.area'],ts);
ar = 0*slp;
ar(1:length(artmp))=artmp;
artmp=[];

xmax = 0.0;  % Maximum streamwise distance found so far

  dx=xyzb(:,1)-x;
  dx=dx.*dx;
  dy=xyzb(:,2)-y;
  dy=dy.*dy;
  dist=sqrt(dx+dy);
  [mindist index]=min(dist);
  %figure(2)
  [x,h,xp,yp,sp,ap]=cfindstrmprofile( xyzb, dir, index, slp, ar );
  if max(x) > xmax, xmax = max(x); end
  plot( x, h )
  hold on

fprintf('Maximum profile distance is %f\n',xmax);








