function [mi,xi,yi]=c2raster(basenm,ts,dx,file_ext)
% C2RASTER: Creates a raster image of a CHILD field, which could be
% elevation or any other field.
%  Usage: [mi xi yi] = c2raster( basenm, ts, npn, {file_ext} )
%  filenm -- name of edge file
%  ts -- time slice to plot
%  npn -- no. of grid points per node
%  file_ext (optional) -- extension of file containing data to be gridded
%  (for example, '.slope' for slope). 
%   G. Tucker, 2008
xyzb=creadxyzb(basenm,ts);
if nargin<4
    m=xyzb(:,3);
else
    m=cread([basenm file_ext],ts);
end
minx=min(xyzb(:,1));
maxx=max(xyzb(:,1));
miny=min(xyzb(:,2));
maxy=max(xyzb(:,2));
[xi yi]=meshgrid((dx/2):dx:maxx,(dx/2):dx:maxy);
m = cinterpclosededges( xyzb(:,1), xyzb(:,2), m, xyzb(:,4) );  % this removes edge effects for plotting
mi=griddata(xyzb(:,1),xyzb(:,2),m,xi,yi);

