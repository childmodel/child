function dd = cfindlayerbase( laydat, nlay, z, laynum )
% CFINDLAYERBASE: Finds the elevation of the bottom of layer NL - LAYNUM at
%                 each node in a CHILD simulation, where NL is the # of layers
%                 at a given node and LAYNUM is the number of layers up from
%                 the bottom we wish to see (0,1,2,...etc). If layer NL is
%                 is always bedrock and the overlying layers sediment, then
%                 using this function with LAYNUM=1 will find the top of
%                 the bedrock-sediment contact.
%  Usage: dd = cfindlayerbase( laydat, nlay, z, laynum )
%   laydat = layer data obtained from READLAYERS
%   nlay = vector of # of layers at each node (also from READLAYERS)
%   z = vector of elevations at each node
%   laynum = # of layers up from bottom we wish to look (0,1,2,...etc)
%       G. Tucker, 1999
nn = size( laydat, 1 );  % number of nodes
dd = z;  % base elevations: start w/ surface elevations, then subtract
nlay = nlay - laynum;  % transform NLAY <- NLAY-LAYNUM
for i=1:nn
  for j=1:nlay(i)  % subtract thickness of each layer up to NL-LAYNUM
    dd(i) = dd(i) - laydat(i,j,1);
  end
end

    
