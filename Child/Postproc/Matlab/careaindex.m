function area_area = careaindex( fname, ts )
% careaindex: computes area accumulation index, the cumulative distribution
% of drainage areas, from a CHILD run of name "fname" and time slice "ts".
%
% Returns an Nx2 matrix in which the first column contains drainage area
% and the second contains the fraction of the total area that has a
% drainage area greater than or equal to the corresponding value in the 
% first column.
% 
% GT, April 2009
%

xyz = creadxyz( fname, ts );
a = cread( [fname '.area'], ts );
n_interior_nodes = length(a);
varea = cread( [fname '.varea'], ts );

area_area = nan( n_interior_nodes, 2 );
area_area(:,1) = a;
area_area(:,2) = varea(1:n_interior_nodes);
area_area = flipud( sortrows( area_area, 1 ) );
area_area(:,2) = cumsum( area_area(:,2) );
area_area(:,2) = area_area(:,2)./max(area_area(:,2));
