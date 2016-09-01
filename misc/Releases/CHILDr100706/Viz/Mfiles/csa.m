function [s,a]=csa( fname, ts, psym )
%CSA: Draws slope-area plot from CHILD simulation.
%   Usage: [s,a] = csa( fname, ts, psym )
%       fname = filename
%       ts = time slice number
%       psym (optional) = plotting symbol, e.g. 'x', 'o', 'v', etc.
%
%  Returns:
%    s: vector of slope values at each interior node
%    a: vector of drainage areas at each interior node
s=cread([fname '.slp'], ts );
a=cread([fname '.area'], ts );
s=s(1:size(a,1));
if nargin==2
  loglog(a/1e6,s,'.')
end
if nargin==3
  loglog(a/1e6,s,psym)
end
set(gca,'fontsize',14);
xlabel( 'Drainage area (km^2)' )
ylabel( 'Slope (m/m)' )
