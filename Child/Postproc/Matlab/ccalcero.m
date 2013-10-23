function erodep=ccalcero( fname, ts, ots )
%CCALCERO: Utility for calculating erosion or sedimentation depth between 
%        time slices in a CHILD run. Assumes that nodes are in the same  
%        order in output files, and that no nodes are added or deleted 
%        during the run. Similar to CCALCERORATE, but with these
%        differences:
%          1. Calculates erosion/deposition by differencing layer
%             thickness rather than surface heights plus uplift rate.
%          2. Calculates thickness change rather than rate.
%
%    Usage: erodep = ccalcero( fname, ts, {ots} )
%      fname = file name without extension
%      ts = time slice, which must be > 1
%      ots = earlier time step (optional; defaults to ts-1)
%
%    Returns: erodep = vector of erosion/sedimentation depth between time  
%             ts and time ots (or ts-1 if ots not specified). To plot a shaded
%             surface colored by erodep, use ctrisurf( fname, ts, erodep )
% 
% GT Aug 2011
if nargin<3
    ots = ts-1;
end

if ots<1, error('ts must be >1'), end
z = cread( [fname '.z'], ts );
[t2,s2,r2] = clayerthick( fname, ts );
[t1,s1,r1] = clayerthick( fname, ots );
erodep = zeros(size(z));
erodep(1:length(t2)) = t2 - t1;



