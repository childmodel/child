function [dist,height]=cfindstrmprofile(xyzb,dir,n)
% CFINDSTRMPROFILE: Traces a stream profile downstream, starting with node n.
%                   Returns the profile in two vectors.
%
%      G. Tucker, 2002
%
dist = [ 0.0 ];
height = [ xyzb(n,3) ];
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
  n = newn;
  if xyzb(n,4) ~= 0, point_in_interior = 0; end
  k=k+1;
  if k>npts, error('Probable endless loop in CFINDSTREAMPROFILE'); end
end