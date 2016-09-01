function [x,h]=cstrmprofile( basenm, ts )
% CSTRMPROFILE: Retrieves and plots a longitudinal stream profile from a
%               CHILD run, starting from a given node which the user selects
%               with the mouse, and moving downstream from there until an
%               open boundary is reached. Returns vectors for the height and
%               distance along the profile.
%
%  Usage: [dist,height] = cstrmprofile( basenm, ts )
%
% Parameters:
%  basenm -- name of files with run
%  ts -- time slice to plot
%
% Functions called: creadxyzb.m, cplotmesh.m, cplotarea.m, cread.m,
%                   cfindstrmprofile.m
%
%   G. Tucker, 2002
%
figure(1)
nd=cplotmesh(basenm,ts);
hold on
ccontour(basenm,ts,4,12);
cplotarea(basenm,nd,ts);
xyzb=creadxyzb(basenm,ts);
dir=cread([basenm '.net'],ts);
dir=dir+1;  % Change of indices from 0..N-1 to 1..N
axis([0 max(xyzb(:,1)) 0 max(xyzb(:,2)) 0 2*max(xyzb(:,3))])
view(0,90)
keepgoing=1;
figure(2)
cla
figure(1)
xmax = 0.0;  % Maximum streamwise distance found so far
fprintf('Use the mouse to select the head of the stream profile.\n')
pt=ginput(1);
fprintf('Point selected: x = %f, y = %f\n',pt(1),pt(2) );
while keepgoing
  dx=xyzb(:,1)-pt(1);
  dx=dx.*dx;
  dy=xyzb(:,2)-pt(2);
  dy=dy.*dy;
  dist=sqrt(dx+dy);
  [mindist index]=min(dist);
  figure(2)
  [x,h]=cfindstrmprofile( xyzb, dir, index );
  if max(x) > xmax, xmax = max(x); end
  plot( x, h )
  hold on
  sz=max(size(x));
  plot( x(sz), h(sz), '.')
  figure(1)
  fprintf('Select another starting point, or press Enter to end.\n');
  pt=ginput(1);
  if isempty(pt)
    keepgoing=0; 
  else
    fprintf('Point selected: x = %f, y = %f\n',pt(1),pt(2) );
  end
end
hold off
figure(2)
hold off
figure(1)
fprintf('Maximum profile distance is %f\n',xmax);








