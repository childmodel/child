function [V,C] = cvoronoi( filenm, ts, c, color_range )
% CVORONOI: Plots a Voronoi diagram from CHILD output.
%
% Usage: [V,C] = cvoronoi( filenm, ts, c, color_range )
%
% GT (prototype), March 2009
%

% Read x,y,z coordinates for desired time step
xyzb = creadxyzb( filenm, ts );

% Handle color variable argument missing
if nargin<3
    c=xyzb(:,3);   % use elevation for color
end

% Handle missing color_range argument and set up color minima and maxima
if nargin<4
    cmax = max(c);
    cmin = min(c);
else
    cmin = color_range(1);
    cmax = color_range(2);
end

% Extract just the x and y coordinates
xy = xyzb(:,1:2);

% Generate the Voronoi diagram data
[V,C] = voronoin(xy);

% Scale the colors
cmap = colormap;
cmin=min(c);
crange = max(c)-cmin;
colortop = size(cmap,1);
color_scale_param = (colortop-1) / crange;
mycolors = round( (c-cmin)*color_scale_param + 1 );
mycolors( mycolors>colortop ) = colortop;
mycolors( mycolors<1 ) = 1;

% For each polygon in the mesh interior, plot the corresponding Voronoi
% cell
num_interior_cells = size( find( xyzb(:,4)<1 ) );
for i=1:num_interior_cells
    
    vertex_indices = C{i};
    %vertex_indices = [vertex_indices vertex_indices(1)]; % wrap around
    
    patch(V(vertex_indices,1),V(vertex_indices,2),mycolors(i),'edgecolor','none');
    %drawnow
    %plot(V(vertex_indices,1),V(vertex_indices,2))
    hold on
    %plot(xy(i,1),xy(i,2),'k.')
    
end

hold off