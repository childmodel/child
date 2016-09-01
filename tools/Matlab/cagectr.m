function di=cagectr( filenm, ts, npn, depth, today )
%CAGECTR: Generates contour plot of ages at given depth from CHILD run.
%  Usage: di = cagectr( filenm, ts, npn, depth, today )
%    filenm -- base name of file to read
%    ts -- time slice number
%    npn -- number of interpolated points per node
%    depth -- depth to look at
%    today -- time (yrs) at end of run, or 0 for max output value
%  G. Tucker, 1998
fprintf( 'Reading layer data...\n' );
lay=readlayers( filenm, ts );
fprintf( 'Reading age data...\n');
dd=cdepthdate( lay, depth );
if today==0, today=max(dd); end
dd=today-dd;
[di xi yi]=tintomatrix( filenm, ts, npn, dd );
di(1,1)=today;
di(1,2)=0;
contourf( xi, yi, di )
colorbar('horiz')
