function [zi,xi,yi]=ccontour(basenm,ts,npn,ncl)
% CCONTOUR: contour plot from Child triangulation
%  Usage: [zi xi yi] = ccontour( basenm, ts, npn, ncl )
%  filenm -- name of edge file
%  ts -- time slice to plot
%  npn -- no. of grid points per node
%  ncl -- no. of contour lines (0 for default)
%   G. Tucker, 1998
xyz=creadxyz(basenm,ts);
minx=min(xyz(:,1));
maxx=max(xyz(:,1));
miny=min(xyz(:,2));
maxy=max(xyz(:,2));
delx=maxx-minx
dely=maxy-miny
dx=sqrt((delx*dely)/(npn*size(xyz,1)));
[xi yi]=meshgrid(0:dx:maxx,0:dx:maxy);
zi=griddata(xyz(:,1),xyz(:,2),xyz(:,3),xi,yi);
%contour(xi,yi,zi,0:5:200,'k')
if ncl==0, contour(xi,yi,zi,'k'), end
if ncl~=0, contour(xi,yi,zi,ncl,'k'),end


