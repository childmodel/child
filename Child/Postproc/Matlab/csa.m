function [s,a]=csa( fname, ts, psym )
%CSA: Draws slope-area plot from CHILD simulation.
%   Usage: [s,a] = csa( fname, ts, psym )
%       fname = filename
%       ts = time slice number
%       psym (optional) = plotting symbol, e.g. 'x', 'o', 'v', etc.
%
s=cread([fname '.slp'], ts );
a=cread([fname '.area'], ts );
if nargin==2
  loglog(a,s(1:size(a,1)),'.')
end
if nargin==3
  loglog(a,s(1:size(a,1)),psym)
end
