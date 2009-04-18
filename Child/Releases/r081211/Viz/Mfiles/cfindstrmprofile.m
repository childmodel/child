function [dist,height,xp,yp,sp,ap]=cfindstrmprofile(xyzb,dir,n,slp,ar)
% CFINDSTRMPROFILE: Traces a stream profile downstream, starting with node n.
%                   Returns the profile in two vectors.
%
%      G. Tucker, 2002; 2008
%
dist = 0.0;
height = xyzb(n,3);
xp = xyzb(n,1);
yp = xyzb(n,2);

if nargin==5
    do_slope_area = 1;
    sp = slp(n);
    ap = ar(n);
else
    do_slope_area = 0;
    sp = 0;
    ap = 0;
end

cumdist = 0.0;
point_in_interior = 1;
k=0;
npts=max(size(xyzb));
while point_in_interior
  newn = dir(n);
  %fprintf('Point %d->%d (%f, %f)\n',n,newn,xyzb(n,1),xyzb(n,2) );
  %pause
  dx = xyzb(newn,1)-xyzb(n,1);
  dy = xyzb(newn,2)-xyzb(n,2);
  edge_length = sqrt( dx*dx + dy*dy );
  cumdist = cumdist + edge_length;
  dist = [ dist cumdist ];
  height = [ height xyzb(newn,3) ];
  xp = [xp xyzb(newn,1)];
  yp = [yp xyzb(newn,2)];
  if do_slope_area
     sp = [sp slp(newn)];
     ap = [ap ar(newn)];
  end
  n = newn;
  if xyzb(n,4) ~= 0, point_in_interior = 0; end
  k=k+1;
  if k>npts, error('Probable endless loop in CFINDSTREAMPROFILE'); end
end