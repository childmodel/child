function di=cexpctr( filenm, ts, npn, depth, maxdate )
%CEXPCTR: Generates contour plot of ages at given depth from CHILD run.
%  Usage: di = cexpctr( filenm, ts, npn, depth, maxdate )
%    filenm -- base name of file to read
%    ts -- time slice number
%    npn -- number of interpolated points per node
%    depth -- depth to look at
%    maxdate -- maximum possible exposure age (use this to calibrate colors)
%  G. Tucker, 1998
fprintf( 'Reading layer data...\n' );
lay=readlayers( filenm, ts );
fprintf( 'Reading exposure age data...\n');
dd=cdepthexp(lay,depth);
[di xi yi]=tintomatrix( filenm, ts, npn, dd );
di(1,1)=0;
di(1,2)=maxdate;
contourf( xi, yi, log10(di+1) )
colorbar('horiz')
