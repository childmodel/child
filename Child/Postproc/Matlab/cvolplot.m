function [t,v] = cvolplot( fname, n, area )
%CVOLPLOT: Plots volume above base level versus time from a Child run. Run must
%          have used OPTTSOUTPUT option.
%    Usage: [t,v] = cvolplot( fname, n, {area} )
%      fname = file name (including path)
%      n = number of lines in <fname>.storm and <fname>.vols files
%      area = total area of simulation mesh (optional)
%    Returns:
%      t = time vector (x-axis)
%      v = volume vector (y-axis)
%   GT, May 2002
stormfile = [fname '.storm'];
volumefile = [fname '.vols'];
sfid = fopen( stormfile, 'r' );
vfid = fopen( volumefile, 'r' );
if sfid <= 0, error(['Unable to open ' stormfile]);end
if vfid <= 0, error(['Unable to open ' volumefile]);end
s = fscanf( sfid, '%f', [3,n]);
t = cumsum( (s(1,:) + s(3,:) ) );   % Row 1 = interstorm dur, row 3 = storm dur
fclose( sfid );
v = fscanf( vfid, '%f', [1,n] );
if nargin==3, v=v./area; end
plot( t, v )
grid on

