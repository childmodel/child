function dd = cdepthdate( laydat, depth )
% CDEPTHDATE: Computes date at a given depth below the surface from
%             a CHILD simulation. We find the layer that encompasses the given
%             depth, and save its date. Returns a vector of dates.
%  Usage: de = cdepthdate( laydat, depth )
%   laydat = layer data obtained from READLAYERS
%   depth = desired depth from which data are extracted
%       G. Tucker, 1998
nn = size( laydat, 1 );  % number of nodes
dd = zeros( nn, 1 );  % date of layer at given depth
for i=1:nn
  cumdep = laydat( i, 1, 1 );  % cumulative depth so far
  lay = 1;
  while cumdep < depth
    lay = lay + 1;
    cumdep = cumdep + laydat( i, lay, 1 );
  end
  dd( i ) = laydat( i, lay, 2 );
end

    
