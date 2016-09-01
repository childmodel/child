function bi=cbasectr( filenm, ts, npn, laynum )
%CBASECTR: Generates contour plot of the elevation at the base of layer
%          N - LAYNUM, where N is the number of layers at a particular point
%          and LAYNUM is a fixed number of layers "up" from the bottom.
%          You can use this function to plot, for example, the bedrock
%          surface beneath a variable number of sediment layers. LAYNUM
%          is coded like:
%              0 - depth to base of lowest layer (N)
%              1 - depth to top of lowest layer (N), eg bedrock surface
%              2 - depth to base of layer N-2
%              3 - etc.
%  Usage: bi = cbasectr( filenm, ts, npn, laynum )
%    filenm -- base name of file to read
%    ts -- time slice number
%    npn -- number of interpolated points per node
%    laynum -- which layer to find base for
%  G. Tucker, Sept. 1999
fprintf( 'Reading layer data...\n' );
[lay nlay]=readlayers( filenm, ts );
fprintf( 'Reading elevation data...\n' );
z=readcnd( [filenm '.z'], ts );
fprintf( 'Finding layer base elevation...\n');
dd=cfindlayerbase( lay, nlay, z, laynum );
ctrisurf(filenm,ts,dd);
fprintf( 'Gridding...\n' );
[bi xi yi]=tintomatrix( filenm, ts, npn, dd );
%contourf( xi, yi, bi )
colorbar('horiz')
