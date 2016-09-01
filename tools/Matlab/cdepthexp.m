function de = cdepthexp( laydat, depth )
% CDEPTHEXP: Computes exposure age at a given depth below the surface from
%            a CHILD simulation. We find the layer that encompasses the given
%            depth, and save its exposure age. Returns a vector of exposure
%            ages.
%  Usage: de = cdepthexp( laydat, depth )
%   laydat = layer data obtained from READLAYERS
%   depth = desired depth from which data are extracted
%              G Tucker, 1998
nn = size( laydat, 1 );  % number of nodes
de = zeros( nn, 1 );  % date of layer at given depth
for i=1:nn
  cumdep = laydat( i, 1, 1 );  % cumulative depth so far
  lay = 1;
  while cumdep < depth
    lay = lay + 1;
    cumdep = cumdep + laydat( i, lay, 1 );
  end
  de( i ) = laydat( i, lay, 3 );
end

    
