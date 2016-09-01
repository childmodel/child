function area_elev = chypsometry( fname, ts )
% chypsometry: computes hypsometric curves, both dimensional and
% normalized, from a CHILD run of name "fname" and time slice "ts".
%
% Returns an Nx4 matrix in which the columns are cumulative area,
% elevation, normalized cumulative area, and normalized elevation.
%
% GT, April 2009
%

xyz = creadxyz( fname, ts );
a = cread( [fname '.area'], ts );  % to get # interior nodes
n_interior_nodes = length(a);
a=[];
varea = cread( [fname '.varea'], ts );
area_elev = nan( n_interior_nodes, 4 );
area_elev(:,1) = varea(1:n_interior_nodes);
area_elev(:,2) = xyz( 1:n_interior_nodes, 3 );
area_elev = flipud( sortrows( area_elev, 2 ) );
area_elev(:,1) = cumsum( area_elev(:,1) );
area_elev(:,3) = area_elev(:,1)./max(area_elev(:,1));
area_elev(:,4) = (area_elev(:,2)-min(area_elev(:,2)))./(max(area_elev(:,2))-min(area_elev(:,2)));
