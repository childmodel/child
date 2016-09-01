function dz=ccalcdz( fname, ts, ots )
%CCALCDZ: Utility for plotting elevation changes between successive
%        time slices in a CHILD run. Assumes that nodes are in the same order 
%        in output files, and that no nodes are added or deleted during the
%        run.
%    Usage: dz = ccalcdz( fname, ts, {ots} )
%      fname = file name without extension
%      ts = time slice, which must be > 1
%      ots = earlier time step (optional; defaults to ts-1)
%    Returns: dz = vector of elevation differences between time ts and
%             time ots (or ts-1 if ots not specified). To plot a shaded
%             surface colored by dz, use ctrisurf( fname, ts, dz )
% 
% GT Dec 2004; modified 5/06
if nargin<3
    ots = ts-1;
end
if ots<1, error('ts must be >1'),end
z1 = cread( [fname '.z'], ts );
z0 = cread( [fname '.z'], ots );
dz = z1-z0;

