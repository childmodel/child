function [t,erorate] = cdvoldtplot( fname, n, area, resamp )
%CDVOLDTPLOT: Plots spatially averaged erosion rates vs time from a Child run,
%          using volume information from the .vols file. Run
%          must have used OPTTSOUTPUT option.
%    Usage: [t,v] = cdvoldtplot( fname, n, area, resamp )
%      fname = file name (including path)
%      n = number of lines in <fname>.storm and <fname>.vols files
%      area (optional) = total area of interior domain, to convert to elev's
%      resamp (optional) = resamples the signal by using only every RESAMP
%                           point (must be an integer)
%    Returns:
%      t = time vector (x-axis)
%      v = - d vol / dt (y-axis)
%   GT, May 2002
stormfile = [fname '.storm'];
volfile = [fname '.vols'];
sfid = fopen( stormfile, 'r' );
vfid = fopen( volfile, 'r' );
if sfid <= 0, error(['Unable to open ' stormfile]);end
if vfid <= 0, error(['Unable to open ' volumefile]);end
s = fscanf( sfid, '%f', [3,n]);
t = cumsum( (s(1,:) + s(3,:) ) );   % Row 1 = interstorm dur, row 3 = storm dur
fclose( sfid );
v = fscanf( vfid, '%f', [1,n] );
fclose( vfid );
if nargin>=3, v=v./area; end
if nargin==4
  t = t(1:resamp:n);
  v = v(1:resamp:n);
end

% Having read the data, now compute derivatives
sz=size(t,2);
dv=diff(v);
if nargin==4
  dt=diff(t);
else
  dt=s(1,n-1)+s(3,n-1);
end
erorate=-dv./dt;
t=t(2:sz);

plot( t, erorate )
grid on

