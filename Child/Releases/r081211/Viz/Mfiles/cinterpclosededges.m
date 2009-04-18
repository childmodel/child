function z = cinterpclosededges( x, y, z, b )
%CINTERPC:PSEDEDGES: CHILD boundary-node elevations are usually fixed at zero.
%This can make for ugly edges in surface plots. This function removes these
%ugly edges by setting the height of each boundary node equal to that of
%its nearest interior neighbor. This is purely cosmetic, of course.
%However, it only does this for "closed" edges.
%
% Usage: z = cinterclosedpedges( x, y, z, b ), where:
%  x, y, z = list of node coordinates (vectors)
%  b = list of node boundary flags (0, 1 or 2) (vector)
% Returns modified z.
%
% GT, August, 2006

% Find the number of nodes and extract x,y coords for interior ones
xi = x( b==0 );  % x-coords of interior nodes
yi = y( b==0 );  % y-coords of interior nodes
closed_nodes = find(b==1);
nn = length(closed_nodes);        % number of nodes

% For each closed boundary node, find nearest interior node and assign z value
for i=1:nn
    dx = xi - x(closed_nodes(i));
    dy = yi - y(closed_nodes(i));
    dr2 = dx.*dx + dy.*dy;  % sq distances to all interior nodes
    [mindr2 minloc] = min(dr2);   % minimum distance and index of node
    z(closed_nodes(i)) = z(minloc);
end