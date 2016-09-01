%CGENOFFSETGRID: Creates an offset triangular grid to be used in a POINTS
%file for a CHILD run.
% Created June 07 GT

% Get starting mesh points
filenm='bmilktest1tm3';
xyzb=creadxyzb(['../bmilktest1/' filenm],1);
x=xyzb(:,1); y=xyzb(:,2); z=xyzb(:,3); b=xyzb(:,4);
outfilenm = 'frankinserttest2.pts';

% Set up bounds and geometry
dx=400*.3048*(1/3);
dy = dx*sind(60);
bottom = 271390+200*.3048;
top = 272860-(400*.3048*(1/3));
%top=bottom;
left = 146180+200*.3048;
right = 147280-(400*.3048*(1/3));
%right = left;

% These are the vectors for x and y points
xn=[];
yn=[];

% Create the high-density points to be added
offsetflag=0;
for yloc=bottom:dy:top

    % Offset every other row by dx/2
    if offsetflag==0
        offsetflag=1;
        offsetamt = 0;
    else
        offsetflag=0;
        offsetamt = dx/2;
    end
    
    % Each point along one row
    for xloc = left:dx:right
        
        xn = [xn xloc+offsetamt];
        yn = [yn yloc];
        
    end
end

% Mask out the original points in the desired area
arearight = max(xn);
arealeft = min(xn);
areatop = max(yn);
areabottom = min(yn);
outer = find( x>arearight | x<arealeft | y>areatop | y<areabottom );
xouter = x(outer); youter = y(outer); zouter = z(outer); bouter = b(outer);

% Interpolate elevations
zn = griddata( x, y, z, xn, yn );

% Combine the low res and high res
n = length(xn);
xall = [rot90(xouter) xn];
yall = [rot90(youter) yn];
zall = [rot90(zouter) zn];
ball = [rot90(bouter) zeros(1,n) ];

% Plot the data for inspection
plot3(xall,yall,zall,'.');

% Write to file
fid = fopen(outfilenm,'w');
fprintf(fid,'%d\n',length(xall));       % Write # of points
for i=1:length(xall)
    fprintf(fid,'%.3f %.3f %.3f %d\n',xall(i),yall(i),zall(i),ball(i) );
end
fclose(fid);


%Notes on what next:
%
% * Read ALL watershed points and extract a desired area
% * Generate and read a GIS mask (irregular shape)