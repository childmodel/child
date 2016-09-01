function [ld,nl]=csetlayerody( fname, ts, kbin, kbout, xmin, ymin, xmax, ymax, numg, format_version )
% CSETLAYERODY: Sets the erodibility value in CHILD layer data. Use this
% when you want to restart a CHILD run using spatially variable
% erodibility. The function reads layer data from a CHILD run, sets the
% erodibility coefficient for all layers within the rectangular region
% between (xmin,ymin) and (xmax,ymax) to kb, and writes the resulting data
% to a file called <fname>_mod.lay<ts-1>
%
% Inputs: fname - base name of file from which to read data
%         ts - time slice to read
%         kbin - erodibility coefficient inside rectangle
%         kbout - erodibility coefficient outside rectangle
%         xmin,ymin - lower left of rectangle
%         xmax,ymax - upper right of rectangle
%         numg (optional) - number of grain-size classes (default=1)
%         format_version (optional) - version of layer format to use (0 or
%                                       1; defaults to 0)
%
% Returns: modified layer data ld, number of layers at each node nl
%
% GT, Nov 2008

% Check numg, default to 1 if not specified
if nargin<9
    numg=1;
end
if nargin<10
    format_version = 0;
end

% Read in initial data (assumes the name of the file is <fname>.lay<ts-1>).
% The layer file will normally have been generated from a previous CHILD
% run.
[ld,nl,time]=creadlayers(fname,ts,numg,format_version);

% Read in the x,y locations for the nodes
xyz = creadxyz(fname,ts);
x = xyz(:,1);
y = xyz(:,2);

% Get a list of index numbers of nodes that fall in the rectangle
nodes_inside_rect = find( x>xmin & x<xmax & y>ymin & y<ymax );

% Set erodibility field
erody = 0*x+kbout;   % This creates erody vector same size as x, all kbout
erody(nodes_inside_rect) = kbin;  % This sets erody to kbin inside the rect

% Write the results to a file
result=cwritelayers(ld, nl, [fname '_mod'], ts, erody, numg, time)



