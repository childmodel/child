function [zi,xi,yi]=ccontour(basenm,ts,npn,ncl)
% CCONTOUR: contour plot from Child triangulation
%  Usage: [zi xi yi] = ccontour( basenm, ts, npn, ncl )
%  filenm -- name of edge file
%  ts -- time slice to plot
%  npn -- no. of grid points per node
%  ncl -- no. of contour lines (0 for default)
%   G. Tucker, 1998; 2008
xyzb=creadxyzb(basenm,ts);
minx=min(xyzb(:,1));
maxx=max(xyzb(:,1));
miny=min(xyzb(:,2));
maxy=max(xyzb(:,2));
delx=maxx-minx;
dely=maxy-miny;
dx=sqrt((delx*dely)/(npn*size(xyzb,1)));
[xi yi]=meshgrid(0:dx:maxx,0:dx:maxy);
z = cinterpedges( xyzb(:,1), xyzb(:,2), xyzb(:,3), xyzb(:,4) );  % this removes edge effects for plotting
zi=griddata(xyzb(:,1),xyzb(:,2),z,xi,yi);
%contour(xi,yi,zi,0:5:200,'k')
if ncl==0, contour(xi,yi,zi,[50:50:400],'k'), end
if ncl~=0, contour(xi,yi,zi,ncl,'k'),end


