function er=ccalcerorate( fname, ts, ots, u )
%CCALCERORATE: Utility for calculating average erosion rates between given
%        time slices in a CHILD run. Assumes that nodes are in the same order 
%        in output files, and that no nodes are added or deleted during the
%        run. Similar to CCALCDZ but with the following differences: it
%        computes and returns the erosion rate, which is equal to:
%
%          E = U - (z(ts)-z(ots))/(t(ts)-t(ots))
%
%        where U is uplift rate over the time period Dt = t(ts) - t(ots).
%        The uplift rate can be a scalar or a vector with the same length
%        as the elevation data.
%
%    Usage: er = ccalcerorate( fname, ts, {ots}, {u} )
%      fname = file name without extension
%      ts = time slice, which must be > 1
%      ots = earlier time step (optional; defaults to ts-1)
%      u = uplift rate(s) (optional; defaults to zero)
%
%    Returns: er = vector of average erosion rates between time ts and
%             time ots (or ts-1 if ots not specified). To plot a shaded
%             surface colored by er, use ctrisurf( fname, ts, er )
% 
% GT Dec 2004; modified 5/06
if nargin<3
    ots = ts-1;
end
if nargin<4
    u=0;
end

if ots<1, error('ts must be >1'), end
[z1,z2,tm1,tm2] = cread2( [fname '.z'], ots, ts );
er = u - (z2-z1)/(tm2-tm1);


