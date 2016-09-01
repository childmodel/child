function [t,v] = cvolplot( fname, n, total_area )
%CVOLPLOT: Plots volume above base level versus time from a Child run. Run must
%          have used OPTTSOUTPUT option.
%    Usage: [t,v] = cvolplot( fname, n )
%      fname = file name (including path)
%      n = number of lines in <fname>.storm and <fname>.vols files
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

% If user gives zero or negative for n, figure it out from file
if nargin<2 || n<=0
   mycomputer = computer;
   if strncmp(mycomputer,'PCWIN',5)
       error('Cannot automatically determine file length on a windows machine.');
   end
   [s,w]=unix(['wc -l ' fname '.vols']);
   n = str2num( w(1:9) );
end

% If user didn't specify area, read it from file
varea = cread([fname '.varea'],1);  % assumes area has stayed constant!
total_area = sum(varea);

s = fscanf( sfid, '%f', [3,n]);
t = cumsum( (s(1,:) + s(3,:) ) );   % Row 1 = interstorm dur, row 3 = storm dur
fclose( sfid );
v = fscanf( vfid, '%f', [1,n] );
fclose( vfid );

% Sometimes the two files won't have equal length, because of output
% buffering or copying of files during an active run. In this case, set the
% length to the shorter of the two (if v is shorter, this is already taken
% care of, so we need only handle the case in which t is shorter than v).
if length(t)<n
    n=length(t);
    v=v(1:n);
end

% Divide volume by total area to get mean altitude and plot
v=v./total_area;
plot( t, v )
grid on

