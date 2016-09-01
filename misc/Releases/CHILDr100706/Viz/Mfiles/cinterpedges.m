function z = cinterpedges( x, y, z, b )
%CINTERPEDGES: CHILD boundary-node elevations are usually fixed at zero.
%This can make for ugly edges in surface plots. This function removes these
%ugly edges by setting the height of each boundary node equal to that of
%its nearest interior neighbor. This is purely cosmetic, of course.
%
% Usage: z = cinterpedges( x, y, z, b ), where:
%  x, y, z = list of node coordinates (vectors)
%  b = list of node boundary flags (0, 1 or 2) (vector)
% Returns modified z.
%
% GT, August, 2006

% Find the number of nodes and extract x,y coords for interior ones
xi = x( b==0 );  % x-coords of interior nodes
yi = y( b==0 );  % y-coords of interior nodes
nn = length(x);        % number of nodes
ni = length(xi);       % number of interior nodes

% For each boundary node, find nearest interior node and assign z value
for i=(ni+1):nn
    dx = xi - x(i);
    dy = yi - y(i);
    dr = sqrt( dx.*dx + dy.*dy );  % distances to all interior nodes
    [mindr minloc] = min(dr);   % minimum distance and index of node
    z(i) = z(minloc);
end