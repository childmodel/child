%arc2child.m
%
% Generates a set of CHILD points from an arcgrid (ascii format) dem.
%
% Francis and Greg, April 2011

% INITIALIZE

% user-controlled parameters
delta = 1;    % child point spacing
demfilename = 'fltr_grd_clp6.txt';
outputfilename = 'frog.pts';
opt_use_llc_coords = false;   % do you want to keep the LLC coords in DEM?

% open boundary method setting (user controlled)
LOWEST_DEM_CELL = 1;
LOWEST_INTERIOR_DEM_CELL = 2;
SPECIFIED_COORDINATES = 3;
open_boundary_method = LOWEST_DEM_CELL;

% read ascii grid DEM
[dem,xllc,yllc,dx,nodataval] = readarc( demfilename );
if ~opt_use_llc_coords
    xllc = 0;
    yllc = 0;
end

% derived parameters
oddrowoffset_x = dx/2 + delta*sind(30);
evenrowoffset_x = dx/2;
bottomrowstart_y = dx/2;
column_spacing = delta;
row_spacing = delta*cosd(30);
sz = size(dem);
nrows = sz(1);
ncols = sz(2);
max_x = xllc + dx*ncols - dx/2;
max_y = yllc + dx*nrows - dx/2;
nhexrows = round( dx*nrows/row_spacing );
nevencols = round( dx*ncols/delta );
noddcols = nevencols - 1;

% total # hex points
if mod(nhexrows,2)==0   % if even #
    nhexpoints = nhexrows*noddcols + (nhexrows/2);
else
    nhexpoints = (nhexrows-1)*noddcols + noddcols + (nhexrows-1)/2;
end

% x,y,z and b
INTERIOR_CELL = 0;
CLOSED_BOUNDARY_CELL = 1;
OPEN_BOUNDARY_CELL = 2;
UNUSED_CELL = 9999;  % code means "this element of the matrix isn't used"
x = zeros(nhexrows,nevencols);
y = zeros(nhexrows,nevencols);
z = zeros(nhexrows,nevencols);
b = UNUSED_CELL+zeros(nhexrows,nevencols);

% "b" defaults to 1 on perimeters
odd_row_of_b = zeros( 1, noddcols );
odd_row_of_b(1) = CLOSED_BOUNDARY_CELL;
odd_row_of_b(end) = CLOSED_BOUNDARY_CELL;
even_row_of_b = zeros( 1, nevencols );
even_row_of_b(1) = CLOSED_BOUNDARY_CELL;
even_row_of_b(end) = CLOSED_BOUNDARY_CELL;

% PROCESS

% Set the x and y coordinates of the hex points
for r = 1:nhexrows
    
   cury = yllc + bottomrowstart_y + (r-1)*row_spacing;
   
   if mod(r,2) == 1   % odd row
       
       newxrow = (xllc + oddrowoffset_x):column_spacing:(max_x-delta/2);
       x(r,1:noddcols) = newxrow;
       newyrow = cury + zeros(size(newxrow));
       y(r,1:noddcols) = newyrow;
       if r>1 && r<nhexrows
           newbrow = odd_row_of_b;
       else
           newbrow = ones( size(odd_row_of_b) );
       end
       b(r,1:noddcols) = newbrow;
       
   else  % even row
       
       newxrow = (xllc + evenrowoffset_x):column_spacing:max_x;
       x(r,:) = newxrow;
       newyrow = cury + zeros(size(newxrow));
       y(r,:) = newyrow;
       
       if r<nhexrows
           newbrow = even_row_of_b;
       else
           newbrow = ones( size(even_row_of_b) );
       end
       b(r,:) = newbrow;
       
   end
   
end

% Set elevations equal to nearest DEM cell, and ...
% Find out their boundary status by checking the nearest cell elevation. If
% it is "no data", set it to 1.
for r=1:nhexrows
    for c=1:nevencols
        
        if b(r,c)~=UNUSED_CELL
            
            demrow = ceil(y(r,c)-yllc)/dx;
            demcol = ceil(x(r,c)-xllc)/dx;
            
            if dem(demrow,demcol) == nodataval
                b(r,c) = CLOSED_BOUNDARY_CELL;
            else
                z(r,c) = dem(demrow,demcol);
            end
            
        end
        
    end
    
end
    
% Get rid of any points flagged as boundaries that have only other
% boundaries as neighbors.
%   First, the 4 corners
if b(2,2)==CLOSED_BOUNDARY_CELL  % lower left
    b(1,1) = UNUSED_CELL;
end
if b(2,nevencols-1)==CLOSED_BOUNDARY_CELL  % lower right
    b(1,noddcols) = UNUSED_CELL;
end
if mod(nhexrows,2)==0   % never need UL or UR when even # rows
    b(nhexrows,1) = UNUSED_CELL;
    b(nhexrows,nevencols) = UNUSED_CELL;
else
    if b(nhexrows-1,2) == CLOSED_BOUNDARY_CELL   % upper left
        b(nhexrows,1) = UNUSED_CELL;
    end
    if b(nhexrows-1,nevencols-1) == CLOSED_BOUNDARY_CELL  % upper right
        b(nhexrows,noddcols) = UNUSED_CELL;
    end
end
%   Second, the bottom row
for c=2:(noddcols-1)
   if b(1,c)==1 && min(b(2,c:(c+1)))~=INTERIOR_CELL
       b(1,c) = UNUSED_CELL;
   end
end
%   Next, top row: the way we handle this depends on whether it's an even
%   or odd row
if mod(nhexrows,2)==0  % if even # rows
   for c=2:(nevencols-1)
       if b(nhexrows,c)==1 && min(b(nhexrows-1,(c-1):c))~=INTERIOR_CELL
        b(nhexrows,c) = UNUSED_CELL;
       end
   end
else
    for c=2:(noddcols-1)
       if b(nhexrows,c)==1 && min(b(nhexrows-1,c:(c+1)))~=INTERIOR_CELL
        b(nhexrows,c) = UNUSED_CELL;
       end        
    end
end
% Now, even rows
for r=2:2:(nhexrows-1)
   % Right side
   if min( b((r-1):(r+1),noddcols) ) ~= INTERIOR_CELL
       b(r,nevencols) = UNUSED_CELL;
   end
   % Left side  
   if min( [ b(r,2) b(r-1,1) b(r+1,1) ] ) ~= INTERIOR_CELL
       b(r,1) = UNUSED_CELL;
   end
   % Interior
   for c=2:(nevencols-1)
      if min( [b(r,c+1) b(r+1,c) b(r+1,c-1) b(r,c-1) b(r-1,c-1) b(r-1,c)] ) ~= INTERIOR_CELL
          b(r,c) = UNUSED_CELL;
      end
   end
end
% And, finally, odd rows
for r=3:2:(nhexrows-1)
   % Right side
   if min( [ b(r,noddcols-1) b(r-1,noddcols) b(r+1,noddcols) ] ) ~= INTERIOR_CELL
       b(r,noddcols) = UNUSED_CELL;
   end
   % Left side  
   if min( [ b(r,2) b(r-1,2) b(r+1,2) ] ) ~= INTERIOR_CELL
       b(r,1) = UNUSED_CELL;
   end
   % Interior
   for c=2:(noddcols-1)
       if min( [b(r,c+1) b(r+1,c+1) b(r+1,c) b(r,c-1) b(r-1,c) b(r-1,c+1)] ) ~= INTERIOR_CELL
           b(r,c) = UNUSED_CELL;
       end
   end
end

% Take the subset of points that are not UNUSED_CELL
xc = x( b~=UNUSED_CELL );
yc = y( b~=UNUSED_CELL );
zc = z( b~=UNUSED_CELL );
bc = b( b~=UNUSED_CELL );

% Assign one or more cells as open boundaries
if open_boundary_method == LOWEST_DEM_CELL
    
    % Find the DEM row and column with the smallest value
    dem(dem<0) = 9999;  % This prevents NODATAVALS from looking like minima
    [minvals,minrows] = min(dem);
    [minval,mincol] = min(minvals);
    minrow = minrows(mincol);
    
    % Find the hex cell closest to this point
    xmin = xllc + dx*mincol - dx/2;
    ymin = yllc + dx*minrow - dx/2;
    dist2 = (xc-xmin).^2 + (yc-ymin).^2;
    [mindist2,minid] = min(dist2);
    
    % Set the boundary code to OPEN_BOUNDARY
    bc(minid) = OPEN_BOUNDARY_CELL;
    
elseif open_boundary_method == LOWEST_INTERIOR_DEM_CELL
    
    % NOT IMPLEMENTED YET
    
else
    
    % NOT IMPLEMENTED YET
    
end

% Set closed boundary elevs to zero (for easier identification in plots; it
% doesn't influence the calculations)
zc( bc==CLOSED_BOUNDARY_CELL ) = 0;


% CLEANUP

% Plot
figure(1)
plot(x(b==INTERIOR_CELL),y(b==INTERIOR_CELL),'b.')
title('Raw hex mesh')
hold on
plot(x(b==CLOSED_BOUNDARY_CELL),y(b==CLOSED_BOUNDARY_CELL),'r.')
plot(x(b==UNUSED_CELL),y(b==UNUSED_CELL),'c.')
legend('Interior cells','Closed boundary cells','Unused cells')
hold off
figure(2)
plot(xc(bc==INTERIOR_CELL),yc(bc==INTERIOR_CELL),'b.')
title('Final hex mesh')
hold on
plot(xc(bc==CLOSED_BOUNDARY_CELL),yc(bc==CLOSED_BOUNDARY_CELL),'r.')
plot(xc(bc==OPEN_BOUNDARY_CELL),yc(bc==OPEN_BOUNDARY_CELL),'g.')
legend('Interior cells','Closed boundary cells','Open boundary cells')
hold off

% Write to file
fid = fopen( outputfilename, 'w' );
fprintf( fid, ' %d\n', length(xc) );
for i=1:length(xc)
   fprintf( fid, '%f %f %f %d\n',xc(i),yc(i),zc(i),bc(i) ); 
end
fclose( fid );
