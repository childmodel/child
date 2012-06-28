function [t,v] = cvolplot( fname, dt )
%CVOLPLOT: Plots volume above base level versus time from a Child run. Run must
%          have used OPTTSOUTPUT option.
%    Usage: [t,v] = cvolplot( fname, n, {dt} )
%      fname = file name
%      n = number of lines in <fname>.storm and <fname>.vols files
%      (optional) dt = time step duration (only if storms aren't random)
%    Returns:
%      t = time vector (x-axis)
%      v = volume vector (y-axis)
%   GT, May 2002

% Open "vols" file
volumefile = [fname '.vols'];
vfid = fopen( volumefile, 'r' );
if vfid <= 0, error(['Unable to open ' volumefile]);end

% If user specified a dt value, assume we don't have random storms
if nargin>2
    use_random_storms = 0;
else
    use_random_storms = 1;
end 

% If applicable, open the "storm" file. If we can't find one, then assume
% we don't have random storms in this run, and use a default dt value
if use_random_storms
    stormfile = [fname '.storm'];
    sfid = fopen( stormfile, 'r' );
    if sfid <= 0
        use_random_storms = 0;
        dt = 1;  % default value
    end
end

% If user gives zero or negative for n, figure it out from file
default_max_number_of_lines = 1e7;
if nargin<2 || n<=0
   mycomputer = computer;
   if strncmp(mycomputer,'PCWIN',5)
       n = default_max_number_of_lines;
   end
   [s,w]=unix(['wc -l ' fname '.vols']);
   n = str2num( w(1:9) );
end

% If user didn't specify area, read it from file
varea = cread([fname '.varea'],1);  % assumes area has stayed constant!
total_area = sum(varea);

% read contents of .vols file
v = fscanf( vfid, '%f', [1,n] );
fclose( vfid );
n = length(v);

% either read contents of storm file, or set "t" to be a series of uniform
% time steps
if use_random_storms
s = fscanf( sfid, '%f', [3,n]);
t = cumsum( (s(1,:) + s(3,:) ) );   % Row 1 = interstorm dur, row 3 = storm dur
fclose( sfid );
else
    t = dt*(1:n);
end

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

