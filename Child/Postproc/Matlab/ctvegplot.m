function [t,v] = ctvegplot( fname, n )
%CVOLPLOT: Plots average fractional veg cover vs time from a Child run. Run
%          must have used OPTTSOUTPUT option.
%    Usage: [t,v] = ctvegplot( fname, n )
%      fname = file name (including path)
%      n = number of lines in <fname>.storm and <fname>.vols files
%    Returns:
%      t = time vector (x-axis)
%      v = fractional vegetation cover vector (y-axis)
%   GT, May 2002
stormfile = [fname '.storm'];
vegfile = [fname '.vcov'];
sfid = fopen( stormfile, 'r' );
vfid = fopen( vegfile, 'r' );
if sfid <= 0, error(['Unable to open ' stormfile]);end
if vfid <= 0, error(['Unable to open ' volumefile]);end
s = fscanf( sfid, '%f', [3,n]);
t = cumsum( (s(1,:) + s(3,:) ) );   % Row 1 = interstorm dur, row 3 = storm dur
fclose( sfid );
v = fscanf( vfid, '%f', [1,n] );
plot( t, v )
grid on

