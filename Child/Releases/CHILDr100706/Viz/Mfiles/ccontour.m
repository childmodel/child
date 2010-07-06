function [zi,xi,yi]=ccontour(basenm,ts,dx,ncl,field,xmin,xmax,ymin,ymax)
% CCONTOUR: contour plot from Child triangulation
%  Usage: [zi xi yi] = ccontour( basenm, ts, npn, ncl )
%  filenm -- name of edge file
%  ts -- time slice to plot
%  npn -- size of grid cells
%  ncl -- no. of contour lines (0 for default)
%   G. Tucker, 1998; 2008
xyzb=creadxyzb(basenm,ts);
if nargin<4
  ncl=0;
end
if nargin<5 || max(field)<=0
  field = xyzb(:,3);
end
if nargin>5
  f=find( xyzb(:,1)>=xmin & xyzb(:,1)<=xmax & xyzb(:,2)>=ymin & xyzb(:,2)<=ymax );
  xyzb = xyzb(f,:);
  field = field(f);
end
minx=min(xyzb(:,1));
maxx=max(xyzb(:,1));
miny=min(xyzb(:,2));
maxy=max(xyzb(:,2));
delx=maxx-minx;
dely=maxy-miny;
%dx=sqrt((delx*dely)/(npn*size(xyzb,1)));
[xi yi]=meshgrid(minx:dx:maxx,miny:dx:maxy);
z = cinterpedges( xyzb(:,1), xyzb(:,2), field, xyzb(:,4) );  % this removes edge effects for plotting
zi=griddata(xyzb(:,1),xyzb(:,2),z,xi,yi);
%contour(xi,yi,zi,0:5:200,'k')
if ncl==0, contour(xi,yi,zi,[50:50:400],'k'), end
if ncl~=0, contour(xi,yi,zi,ncl,'k'),end


