function [outlet_distance,node_count,dist] = cwidthfn( fname, ts, nbins, maxlen )
% cwidthfn: computes the width function from a CHILD run of name 
% "fname" and time slice "ts". The optional "nbins" parameter allows you to
% specify the number of bins for the frequency distribution. The maxlen
% parameter allows you to set the maximum length in for the histogram.
%
% Returns:
%  outlet_distance: vector with outlet distances for each interior node
%  node_count: number of nodes within a given distance increment of outlet
%  dist: distance increments (bin centers) in histogram
%
% GT, April 2009
%

% Read 'net' file containing id #s of downstream neighbor nodes
downstrm_nbrs = cread( [fname '.net'], ts );
downstrm_nbrs = downstrm_nbrs+1;  % shift from 0..N-1 to 1..N numbering
n_interior_nodes = length(downstrm_nbrs);

% Read x,y,z and boundary info
xyzb = creadxyzb( fname, ts );
boundary = xyzb(:,4);

% Initialize matrices
outlet_distance = nan( size(downstrm_nbrs) );

% For each node, step downstream from node to node until an outlet is
% reached, calculating the total distance along the way.
for i=1:n_interior_nodes
   k=i;
   endless_loop_check = 0;
   dist_so_far = 0;
   x = xyzb(i,1);
   y = xyzb(i,2);
   
   % Keep going as long as we're in the interior
   while boundary(k)==0
       downstrm_id = downstrm_nbrs(k);
       x_down = xyzb( downstrm_id, 1 );
       y_down = xyzb( downstrm_id, 2 );
       dist_so_far = dist_so_far + sqrt( (x-x_down)^2 + (y-y_down)^2 );
       k = downstrm_id;
       x = x_down;
       y = y_down;
       endless_loop_check = endless_loop_check+1;
       if endless_loop_check>n_interior_nodes
           error('Too many iterations in WHILE loop');
       end
   end
   
   outlet_distance(i) = dist_so_far;
   
end

if nargin>2 && nbins>0
    if nargin>3 && maxlen>0
        binwid = maxlen/nbins;
        [node_count,dist] = hist(outlet_distance,0.5*binwid:binwid:maxlen);
    else
        [node_count,dist] = hist(outlet_distance,nbins);
    end
else
    avg_area = mean(varea(1:n_interior_nodes));
    char_length = sqrt( avg_area );
    [node_count,dist] = hist(outlet_distance,...
        0.5*char_length:char_length:max(outlet_distance));
end



