function zi=ccontourxyc(npn,ncl,xyz,c)
% CCONTOUR: contour plot from Child triangulation
%  npn -- no. of grid points per node
%  ncl -- no. of contour lines (0 for default)
%  xyz -- xyz matrix
%  c -- field to plot (0 for elevation)
if c==0, c=xyz(:,3); end
minx=min(xyz(:,1));
maxx=max(xyz(:,1));
miny=min(xyz(:,2));
maxy=max(xyz(:,2));
delx=maxx-minx
dely=maxy-miny
dx=sqrt((delx*dely)/(npn*size(xyz,1)))
[xi yi]=meshgrid(0:dx:maxx,0:dx:maxy);
zi=griddata(xyz(:,1),xyz(:,2),c,xi,yi);
if ncl==0, contour(xi,yi,zi), end
if ncl~=0, contour(xi,yi,zi,ncl),end


